#include <iostream>
#include <map>

#include <sam.h>

int main(int argc, char ** argv)
{
    if (argc != 2) {
        std::cerr << "Usage: pac_bam_scan FILE.bam\n";
        return 1;
    }

    samFile * samfile = sam_open(argv[1], "r");
    bam_hdr_t * hdr = sam_hdr_read(samfile);

    bam_hdr_destroy(hdr);
    sam_close(samfile);

    return 0;
}
