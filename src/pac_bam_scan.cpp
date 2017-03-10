#include <algorithm>
#include <cstdio>
#include <iostream>
#include <map>
#include <vector>
#include <clocale>

#include <sam.h>



struct ChunkStats
{
    ChunkStats() = default;

    ChunkStats(int length, int num_mismatches, int num_insertions, int num_deletions) :
        length(length), num_mismatches(num_mismatches), num_insertions(num_insertions),
        num_deletions(num_deletions)
    {}

    int length { 0 };
    int num_mismatches { 0 };
    int num_insertions { 0 };
    int num_deletions { 0 };
};


struct ReadAlignmentStats
{
    int readLength { 0 };
    std::vector<ChunkStats> chunkStats;
    std::vector<std::pair<int, int> > alignedRanges;

    void registerAlignment(bam1_t * record);

    double mmRate() const {
        int lenSum = 0;
        int mmSum = 0;
        for (auto const & stats : chunkStats) {
            lenSum += stats.length;
            mmSum += stats.num_mismatches;
        }
        if (lenSum == 0)
            return 0;

        return ((double)mmSum) / lenSum * 100.0;
    }

    double insRate() const {
        int lenSum = 0;
        int insSum = 0;
        for (auto const & stats : chunkStats) {
            lenSum += stats.length;
            insSum += stats.num_insertions;
        }
        if (lenSum == 0)
            return 0;

        return ((double)insSum) / lenSum * 100.0;
    }

    double delRate() const {
        int lenSum = 0;
        int delSum = 0;
        for (auto const & stats : chunkStats) {
            lenSum += stats.length;
            delSum += stats.num_deletions;
        }
        if (lenSum == 0)
            return 0;

        return ((double)delSum) / lenSum * 100.0;
    }

    int alignedLength() const
    {
        std::vector<std::pair<int, int> > sortedRanges(alignedRanges);
        std::sort(sortedRanges.begin(), sortedRanges.end());

        int current = 0;
        int sum = 0;
        for (auto const & range : sortedRanges) {
            int first = std::max(range.first, current);
            if (range.second > first)
                sum += range.second - first;
            current = range.second;
        }
        return sum;
    }

    double alignedPercentage() const
    {
        if (alignedRanges.empty())
            return 0;
        else
            return (alignedLength() / static_cast<double>(readLength)) * 100;
    }

private:

    ChunkStats computeChunkStats(bam1_t * record) const
    {
        int numMismatches = countMismatches(record);

        int length = 0;
        int numInsertions = 0;
        int numDeletions = 0;

        int start = -1;
        int end = -1;
        int current = 0;

        uint32_t * cigar = bam_get_cigar(record);
        if (!(record->core.flag & BAM_FUNMAP))
            for (int i = 0; i < record->core.n_cigar; ++i) {
                uint32_t cigarElem = cigar[i];

                switch (bam_cigar_opchr(cigarElem)) {
                    // skipping of read bases
                    case 'S':
                    case 'H':
                        current += bam_cigar_oplen(cigarElem);
                    break;

                    // base in ref but not read
                    case 'D':
                        numDeletions += bam_cigar_oplen(cigarElem);
                        break;

                    // "using up" read bases
                    case 'I':
                        numInsertions += bam_cigar_oplen(cigarElem);
                    case 'M':
                    case '=':
                    case 'X':
                        if (start == -1)
                            start = current;
                        current += bam_cigar_oplen(cigarElem);
                        end = current;
                    break;
                }
            }

        return ChunkStats(end - start, numMismatches, numInsertions, numDeletions);
    }

    int countMismatches(bam1_t * record) const
    {
        if (record->core.flag & BAM_FUNMAP)
            return 0;

        int numMismatches = 0;

        uint8_t * md = bam_aux_get(record, "MD");
        if (md == NULL)
            return 0;
        if (md[0] != 'Z')
            throw std::runtime_error("Unexpected MD");
        else
            ++md;

        enum { MATCH, MISMATCH, DEL } state = MATCH;

        int start = 0;
        int current = 0;
        uint8_t * ptr = md;
        for (; *ptr != '\0'; ++ptr) {
            if (isdigit((char)*ptr)) {
                if (state == MATCH) {
                    current += 1;
                } else {
                    // ignore DELs, we get this from CIGAR, only handle MISMATCH
                    if (state == MISMATCH) {
                        numMismatches += (current - start);
                    }

                    // begin MATCH state with next
                    state = MATCH;
                    start = current;
                    current += 1;
                }
            } else if (*ptr == '^') {
                // ignore MATCH, we count mismatches only
                if (state == MISMATCH) {
                    numMismatches += (current - start);
                }
                // register DEL state but don't handle
                state = DEL;
                start = current;
                current += 1;
            } else {
                // letter, only need to switch of in case of MATCH
                if (state == MATCH) {
                    state = MISMATCH;
                    start = current;
                }
                current += 1;
            }
        }
        if (state == MISMATCH)
            numMismatches += (current - start);

        return numMismatches;
    }

    std::pair<int, int> alignedRange(bam1_t * record) const
    {
        if (record->core.flag & BAM_FUNMAP)
            return std::make_pair(0, -1);  // unmapped

        int start = -1;
        int end = -1;
        int current = 0;

        uint32_t * cigar = bam_get_cigar(record);
        for (int i = 0; i < record->core.n_cigar; ++i) {
            uint32_t cigarElem = cigar[i];

            switch (bam_cigar_opchr(cigarElem)) {
                // skipping of read bases
                case 'S':
                case 'H':
                    current += bam_cigar_oplen(cigarElem);
                break;

                // "using up" read bases
                case 'M':
                case 'I':
                case '=':
                case 'X':
                    if (start == -1)
                        start = current;
                    current += bam_cigar_oplen(cigarElem);
                    end = current;
                break;
            }
        }

        return std::make_pair(start, end);
    }
};


void ReadAlignmentStats::registerAlignment(bam1_t * record)
{
    // Update read length
    readLength = std::max(readLength, record->core.l_qseq);

    if (record->core.flag & BAM_FUNMAP)
        return;

    // Get chunk length
    alignedRanges.push_back(alignedRange(record));

    // Compute number of insertions, deletions, and mismatches.
    chunkStats.push_back(computeChunkStats(record));
}


class AlignmentProcessor
{
public:
    AlignmentProcessor(bam_hdr_t * hdr) : hdr(hdr)
    {}

    void registerAlignment(bam1_t * record);

    void print(std::ostream & out) const;

private:

    // pointer to the BAM header
    bam_hdr_t * hdr { 0 };

    // statistics for each read
    std::map<std::string, ReadAlignmentStats> readStats;
};


void AlignmentProcessor::registerAlignment(bam1_t * record)
{
    std::string qname(bam_get_qname(record));
    auto & stats = readStats[qname];
    stats.registerAlignment(record);
}

void AlignmentProcessor::print(std::ostream & out) const
{
    out << readStats.size() << " reads with stats\n";

    std::map<int, int> lengthHistogram;
    for (auto const p : readStats)
        lengthHistogram[p.second.readLength] += 1;
    for (auto const & h : lengthHistogram)
        out << "LEN\t" << h.first << "\t" << h.second << "\n";

    std::vector<char> buffer(1000);

    out << "HREAD\tname\treadLength\taliLength\taliPerc\tmmRate\tinsRate\tdelRate\n";
    for (auto const & p : readStats) {
        snprintf(&buffer[0], 9999, "READ\t%s\t%d\t%d\t%6.2f\t%6.2f\t%6.2f\t%6.2f\n",
                 p.first.c_str(), p.second.readLength, p.second.alignedLength(),
                 p.second.alignedPercentage(),
                 p.second.mmRate(), p.second.insRate(), p.second.delRate());
        out << &buffer[0];
    }
}


int main(int argc, char ** argv)
{
    if (argc != 2) {
        std::cerr << "Usage: pac_bam_scan FILE.bam\n";
        return 1;
    }

    std::cerr << "Opening " << argv[1] << "\n";
    samFile * samfile = sam_open(argv[1], "r");
    bam_hdr_t * hdr = sam_hdr_read(samfile);
    bam1_t * rec = bam_init1();

    AlignmentProcessor proc(hdr);

    for (uint64_t no = 0; sam_read1(samfile, hdr, rec) >= 0; ++no)
    {
        // TODO: ignore non 1..22, X, Y alignment
        if (no % 10000 == 0) {
            if (rec->core.tid >= 0) {
                setlocale(LC_NUMERIC, "");
                std::vector<char> buffer(100);
                snprintf(&buffer[0], 99, "%'d", rec->core.pos + 1);
                std::cerr << "currently at " << hdr->target_name[rec->core.tid] << ":" << &buffer[0] << "\n";
            } else {
                std::cerr << "currently at unaligned read\n";
            }
        }

        proc.registerAlignment(rec);
    }

    proc.print(std::cout);

    std::cerr << "Closing file again\n";
    bam_destroy1(rec);
    bam_hdr_destroy(hdr);
    sam_close(samfile);

    return 0;
}
