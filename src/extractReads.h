#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include <iostream>
#include <seqan/align.h>
#include <seqan/bam_io.h>
#include <iostream>
#include <fstream>
#include <seqan/find.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <map>
#include <iostream>
#include <fstream>


using namespace seqan;
using namespace std;



std::vector<BamAlignmentRecord > extractReads(Options &o)
{
    std::vector<BamAlignmentRecord> uniqueExtractedReads;
    std::cout << "Open bam file\n";
    BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(o.readsFile)))
    {
        std::cerr << "ERROR: Could not open " << o.readsFile << std::endl;
        exit(1);
    }
    BamHeader header;
    std::vector<std::vector<BamAlignmentRecord > > recordtable;
    try
    {
        // Copy header.
        readHeader(header, bamFileIn);
        // Copy records.
        BamAlignmentRecord record;
        while (!atEnd(bamFileIn))
        {
            readRecord(record, bamFileIn);
            if(hasFlagUnmapped(record))
                continue;
            if(hasFlagSecondary(record))
                continue;
            uniqueExtractedReads.push_back(std::move(record));
        }
    }
    catch (Exception const & e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        exit(1);
    }
    close(bamFileIn);

    std::cout << "Loaded " << uniqueExtractedReads.size() << " reads.\n";

    return uniqueExtractedReads;

}
