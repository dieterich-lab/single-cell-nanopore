/*==================================================

   Flexbar - flexible barcode and adapter removal

   Version 3.5.0

   BSD 3-Clause License

   uses SeqAn library release 2.4.0
   and TBB library 4.0 or later


             Developer:  Johannes Roehr

   Former contributors:  Matthias Dodt
                         Benjamin Menkuec
                         Sebastian Roskosch


   https://github.com/seqan/flexbar

===================================================*/

#include "Flexbar.h"
#include "extractReads.h"
#include "removeCDNA.h"
#include "FlexbarTypes.h"
#include <chrono>

int main(int argc, const char* argv[]){

    //TODO set Log file to leftTail in the beginning
	using namespace std;
	using namespace seqan;

    auto start = std::chrono::high_resolution_clock::now();

	const string version = "3.5.0";
	const string date    = "May 2019";

	ArgumentParser parser("singleCellPipe");

	defineOptions(parser, version, date);
	parseCmdLine(parser, version, argc, argv);

	Options o;
	initOptions(o, parser);
    loadOptions(o, parser);

    AlignmentResults res;

    {
        std::cout << "Extract Reads\n";
        std::vector<BamAlignmentRecord > recordstable = extractReads(o);
//             std::cout << recordstable.size() << " finished\n";
        std::cout << "Remove cDNA from Reads\n";
        removeCDNA(o,recordstable);
    }

    //TODO use output prefix // maybe target
    o.fstrmOut.close();

    std::string homopolymersRighttmp = o.htrimRight;
//         std::string homopolymersLefttmp = rc_string(homopolymersRighttmp);

    o.htrimRight = "";

    std::string targetName = o.targetName;
//     o.targetName += "_P1";
    std::string logFileName = o.targetName + ".res";
    openOutputFile(o.fstrmOut, logFileName);
    o.out = &o.fstrmOut;
    //TODO //overwrite adapter parameters in options
    //TODO use own parameters //overwrite adapter parameters in options

    //P1 alignment

    //use adapter para
//     *o.out << "Primer Alignment\n";
    std::cout << "Primer Alignment\n";
    {
        int tmpbundleSize = o.bundleSize;
        o.bundleSize = 256;
        startComputation(o, res);
        o.bundleSize = tmpbundleSize;
    }

    std::cout << "Finished\n";

    o.skipOutput = true;
    o.barcodeAlignment = true;
    o.tableHeader = false;

/*
    typedef std::map<seqan::CharString, uint32_t>  PScore;
    std::cout << "Print leftTail: \n";
    for(int i = 0; i < res.leftTail.size(); ++i){
        std::cout << get<0>(res.leftTail[i]) << "\n";
        std::cout << "map size: " << res.leftTailScores.size() << "\n";
        PScore::iterator it = res.leftTailScores.find(get<0>(res.leftTail[i]));
        if(it != res.leftTailScores.end())
        {
            std::cout << it->second << "\n";
        }
        std::cout << get<1>(res.leftTail[i]) << "\t";
        std::cout << get<2>(res.leftTail[i]) << "\n";
    }*/


    appendOutputFile(o.fstrmOut, logFileName);
    o.out = &o.fstrmOut;
//     *o.out << "Left Barcode Alignment\n";

    o.adapterSeq = "";
    o.useAdapterFile = true;
    o.adapterFile = o.whitelist;
//         std::cout << "Adapter file: " << o.adapterFile << "\n";

    o.fastaRecords = std::move(res.leftTail);

        //TODO always check if it is empty

    o.logEverything = true;
//     o.skipStatistics = false;
    o.a_match = o.barcode_match;
    o.a_mismatch = o.barcode_mismatch;
    o.a_gapCost = o.barcode_gapCost;
    o.a_errorRate = o.barcode_errorRate;
    o.a_end = flexbar::TrimEnd::LTAIL;
    o.rcMode    = flexbar::RevCompMode::RCOFF;
    o.a_min_overlap = o.b_min_overlap;
    o.htrimRight = homopolymersRighttmp;

    std::cout << "Left Tail Size: " << o.fastaRecords.size() << "\n";
    {
        startComputation(o, res);
    }
    std::cout << "Finished Left Tail Barcode Alignment\n";

    o.fastaRecords = std::move(res.rightTail);
    std::cout << "Right Tail Size: " << o.fastaRecords.size() << "\n";
    appendOutputFile(o.fstrmOut, logFileName);
    o.out = &o.fstrmOut;

//     *o.out << "Right Barcode Alignment\n";

    //comp
    o.a_end = flexbar::TrimEnd::RTAIL;
    o.rcMode = flexbar::RevCompMode::RCONLY;

    {
        startComputation(o, res);
    }
    std::cout << "Finished Right Tail Barcode Alignment\n";

    auto end = std::chrono::high_resolution_clock::now();
//         std::chrono::duration<int, std::kilo> elapsed = end - start;
    std::cout << "Total wall clock runtime " <<chrono::duration_cast<chrono::seconds>(end-start).count() << "s." << "\n";








	return 0;
}
