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
#include "splitReads.h"
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

        //TODO extract Reads
        //additionally make sure they are unique
        //TODO std::move??
        {
            std::cout << "Extract Reads\n";
            std::vector<BamAlignmentRecord > recordstable = extractReads(o);
//             std::cout << recordstable.size() << " finished\n";
            std::cout << "Remove cDNA from Reads\n";
            removeCDNA(o,recordstable);
        }

//         for(int i = 0; i < 3; ++i){
//             std::cout << o.fastaRecords[i].first << "\t" << o.fastaRecords[i].second << "\n";
//         }

//         o.fastaRecords.erase(o.fastaRecords.begin() + 15, o.fastaRecords.end());
        //TODO use output prefix // maybe target
        o.fstrmOut.close();

        std::string targetName = o.targetName;
        o.targetName += "_P1";
        std::string logFileName = o.targetName + ".log";
        openOutputFile(o.fstrmOut, logFileName);
        o.out = &o.fstrmOut;
        //TODO //overwrite adapter parameters in options
        //TODO use own parameters //overwrite adapter parameters in options

        //P1 alignment

        //use adapter para
        std::cout << "Primer Alignment\n";
        {
            int tmpbundleSize = o.bundleSize;
            o.bundleSize = 256;
            startComputation(o);
            o.bundleSize = tmpbundleSize;
        }

    std::cout << "Finished\n";

    o.skipOutput = true;

    splitReads(o);
/*
        std::cout << "left: " << o.leftTail.size() << "\n\n";
        for(int i = 0; i < 3; ++i){
            std::cout << o.leftTail[i].first << "\t" << o.leftTail[i].second << "\n";
        }

        std::cout << "right: " << o.rightTail.size() << "\n\n";
        for(int i = 0; i < 3; ++i){
            std::cout << o.rightTail[i].first << "\t" << o.rightTail[i].second << "\n";
        }*/

        o.fstrmOut.close();
//         exit(0);
        o.targetName = targetName + "_leftTail";
        logFileName = o.targetName + ".log";
        appendOutputFile(o.fstrmOut, logFileName);
        o.out = &o.fstrmOut;

        o.adapterSeq = "";
        o.useAdapterFile = true;
        o.adapterFile = o.whitelist;
//         std::cout << "Adapter file: " << o.adapterFile << "\n";

        o.fastaRecords = std::move(o.leftTail);

        //TODO always check if it is empty


//         if(o.bundleSize == 256)
        o.bundleSize = 4;
        o.logEverything = true;
        o.a_match = o.barcode_match;
        o.a_mismatch = o.barcode_mismatch;
        o.a_gapCost = o.barcode_gapCost;
        o.a_errorRate = o.barcode_errorRate;
        o.a_end = flexbar::TrimEnd::LTAIL;
        o.rcMode    = flexbar::RevCompMode::RCOFF;
        o.a_min_overlap = o.b_min_overlap;

        {
            startComputation(o);
        }
        std::cout << "Finished Left Tail Barcode Alignment: " << logFileName << "\n";

        o.fastaRecords = std::move(o.rightTail);
        o.targetName = targetName + "_rightTail";
        logFileName = o.targetName + ".log";
        openOutputFile(o.fstrmOut, logFileName);
        o.out = &o.fstrmOut;

        //comp
        o.a_end = flexbar::TrimEnd::RTAIL;
        o.rcMode = flexbar::RevCompMode::RCONLY;

        {
            startComputation(o);
        }
        std::cout << "Finished Right Tail Barcode Alignment: " << logFileName << "\n";

        auto end = std::chrono::high_resolution_clock::now();
//         std::chrono::duration<int, std::kilo> elapsed = end - start;
        std::cout << "Total wall clock runtime " <<chrono::duration_cast<chrono::seconds>(end-start).count() << "s." << "\n";








	return 0;
}
