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
    bool logregions = true;
    if(o.regionsLog.compare("") == 0)
        logregions = false;
    //TODO add as parameters
    int threshold = -1;
    int step = 1;
    int overlap = 0;
    bool verbose = false;
//     bool overlapping = false;
    bool rmDup = true;
    int threads = o.nThreads;
    int mergeRegionsDis = 1000;


    //Read gtf file
    std::cout << "Open regions file\n";
    ifstream inputFile(toCString(o.regionsFile));
    string line;
    vector<std::tuple<CharString, uint32_t, uint32_t, CharString> > table/*completeTable*/;
    vector<string> oldcontigs;
    if(inputFile.is_open()){
//         bool head = true;
        while (getline(inputFile, line, '\n'))
        {
            if(line.length() < 10)
                continue;
//             if(head){
//                 head = false;
//                 continue;
//             }
            std::istringstream sline(line);
            std::string rowelement;
            std::vector<string> tmprow;
            int k = 0;
            while (getline(sline, rowelement, '\t')){
//                 std::cout << "k: " << k << "\t" << rowelement << "\n";
                if(k == 0 || k == 3 || k == 4)
                    tmprow.push_back(rowelement);
                if(k == 8){
                    std::istringstream sname(rowelement);
                    std::string genename;
                    if(getline(sname, genename, ';'))
                        tmprow.push_back(genename);
                    break;
                }
                ++k;
            }
            CharString refid = tmprow[0];
            // -1 to convert from bam start coordinate 1 to 0 coordinate system in seqan
            uint32_t s = std::stoi(tmprow[1]) - overlap - 1;
            uint32_t e = std::stoi(tmprow[2]) + overlap - 1;
            s = (s < 0) ? 0 : s;
            CharString name = tmprow[3];
            table.push_back(std::make_tuple(refid, s, e, name)); // completeTable
        }
//         std::cout << "Finished reading\n";

//         for (unsigned i = 0; i < table.size(); i++){
//             std::cout << std::get<0>(table[i]) << "\t" << std::get<1>(table[i]) << "\t" << std::get<2>(table[i]) << "\n";
//         }
//         std::cout << "\n";
        inputFile.close();
    }
    else
    {
        std::cout << "Could not open file: " << o.regionsFile << "\n";
    }


    // merge close regions
    for(int i = 0; i < table.size() - 1; ++i){
        if(get<0>(table[i]) == get<0>(table[i + 1])){
            if(get<2>(table[i]) > get<1>(table[i + 1]) - mergeRegionsDis)
            {
                get<2>(table[i]) = (get<2>(table[i]) > get<2>(table[i + 1])) ? get<2>(table[i]) : get<2>(table[i + 1]);
                get<3>(table[i]) += ", ";
                get<3>(table[i]) += get<3>(table[i + 1]);
                table.erase(table.begin() + i + 1);
                --i;
            }
        }
    }

    // check contigs are sorted in some order
    string lastContig = "";
    vector<string> scannedContigs;
    for(int i = 0; i < table.size(); ++i)
    {
        if(lastContig.compare(toCString(get<0>(table[i]))) != 0){
            lastContig = toCString(get<0>(table[i]));

            for(int j = 0; j < scannedContigs.size() && !scannedContigs.empty(); ++j){
//                 std::cout << "l: " << lastContig << "\tscan\t" << scannedContigs[j] << "\n";
                if(lastContig.compare(scannedContigs[j]) == 0){
                    std::cerr << "Contigs in the gtf file are not sorted\n";
                    exit(0);
                }
            }
            scannedContigs.push_back(lastContig);
        }
    }


//     for (unsigned i = 0; i < table.size(); i++){
//         std::cout << std::get<0>(table[i]) << "\t" << std::get<1>(table[i]) << "\t" << std::get<2>(table[i]) << "\n";
//     }
//     std::cout << "\n";


    // Open input file, BamFileIn can read SAM and BAM files.

 /*
    if(completeTable.size() < interval)
        interval = completeTable.size();



    int p = 0;
    while(completeTable.size() > 0){
        vector<std::tuple<CharString, uint32_t, uint32_t> > table(interval);
        if(completeTable.size() > interval){
            table.insert(table.begin(), completeTable.begin(), completeTable.begin() + interval);
            completeTable.erase(completeTable.begin(), completeTable.begin() + interval);
        }
        else
        {
            table = completeTable;
            completeTable.clear();
        }*/

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
            recordtable.resize(table.size());

            string lastContig = "ou";
            int st = 0;
            int end = 0;
            bool stNew = true;
            bool endNew = false;
            bool skip = false;
            uint32_t k = 0;
            while (!atEnd(bamFileIn))
            {
                readRecord(record, bamFileIn);
                //use map to jump to correct chromosom //use start pos and length
                if(hasFlagUnmapped(record))
                    continue;
                if(hasFlagSecondary(record) && !o.rmMulti)
                    continue;
                string recordContig = toCString(getContigName(record, bamFileIn));
                uint32_t recordBegin = record.beginPos;
                uint32_t recordEnd = recordBegin + getAlignmentLengthInRef(record); //length(record.seq);
                //determine search range
    //             bool same = true;
                if(lastContig.compare(recordContig) != 0){
                    if(verbose)
                        std::cout << "Calc new Range for: \n";
                    if(!skip)
                    {
//                         table.erase(table.begin() + st, table.begin() + end);
                        for(int i = st; i < end; ++i)
                            table.erase(table.begin() + st);

                    }
    //                 same = false;
                    skip = false;
                    stNew = true;
                    string rowContig;
                    for(int i = 0; i < table.size(); ++i){
                        rowContig = toCString(std::get<0>(table[i]));
                        if(recordContig.compare(rowContig) == 0 && stNew){
                            stNew = false;
                            st = i;
                        }
                        if(!stNew && recordContig.compare(rowContig) != 0)
                        {
                            end = i;
                            break;
                        }
                    }

                    // in case last element matches
                    if(recordContig.compare(rowContig) == 0)
                        end = table.size();

                    lastContig = recordContig;
                    if(verbose){
                        std::cout << "Chrom: " << lastContig << "\n";
                        std::cout << "Start: " << st << "\tEnd: " << end << "\n";
                        std::cout << "Record: " << k << "\n";
                    }

                    // in case no element matches
                    if(stNew){
                        std::cout << "Skip not in gtf file: " << lastContig << "\n";

//                         string finishedcontig = toCString(std::get<0>(table[st]));
                        for(int i = 0; i < oldcontigs.size(); ++i){
//                             std::cout << "Check: " << oldcontigs[i] << "\told: " << lastContig << "\n";
                            if(lastContig.compare(oldcontigs[i]) == 0)
                            {
                                //TODO add option to ignore unsorted bamfile dont erase then
                                std::cerr << "Error records in Bamfile are not sorted\n"; // read extraction will be slow
                                exit(0);
                            }
                        }


                        skip = true;
                    }
                    oldcontigs.push_back(lastContig);

                }
                ++k;

                if(skip){
                    st = table.size();
                    end = table.size();
                    continue;
                }

                //TODO assume sorted
                #pragma omp parallel for num_threads(threads) schedule(static)
                for(int i = st; i < end; ++i){
                    //check row
                    string rowContig = toCString(std::get<0>(table[i]));
                    uint32_t rowBegin = std::get<1>(table[i]);
                    uint32_t rowEnd = std::get<2>(table[i]);

//                     std::cout << "rBegin: " << recordBegin  << "\t" << rowBegin << "\t" << recordEnd << "\n";
                    // check in nanopore read ends or starts inside annotated region
                    if((recordBegin >= rowBegin && recordBegin <= rowEnd) || (recordEnd >= rowBegin && recordEnd <= rowEnd))
                    {
//                         std::cout << "accepted" << "\n";
                        recordtable[i].push_back(record);
                    }
                    /*
                    // check if nanopore reads is overlapping
                    else if (overlapping && recordBegin <= rowBegin && recordEnd >= rowEnd)
                    {
                        if(threshold == -1){
                                recordtable[i].push_back(record);
                        }
                        else if((recordBegin + threshold > rowBegin) && (recordEnd <= rowEnd + threshold))
                        {
                                recordtable[i].push_back(record);
                        }
                    }*/
                }
            }


        }
        catch (Exception const & e)
        {
            std::cerr << "ERROR: " << e.what() << std::endl;
            exit(1);
        }
        close(bamFileIn);
//     }


    typedef map<seqan::CharString, short>  IDMAP;
    IDMAP idMap_occ;

    if (o.rmMulti){
        for(int i = 0; i < recordtable.size(); ++i){
            for(int j = 0; j < recordtable[i].size(); ++j){
                auto & record = recordtable[i][j];
                IDMAP::iterator it;
                it = idMap_occ.find(record.qName);
                if(it != idMap_occ.end())
                {
//                     std::cout << (int) it->second << "\n";
                    if (it->second == 1)
                        idMap_occ[record.qName] = 2;
                }
                else
                {
                    idMap_occ[record.qName] = 1;
                }
            }
        }
    }

    std::vector<BamAlignmentRecord> uniqueExtractedReads;
    //post processing making reads unique
    map<seqan::CharString, short> idMap;
    uint32_t dups = 0;
    uint32_t empty_dups = 0;
    uint32_t unique = 0;

    ofstream regionsLog;
    if(logregions || true)
    {
    regionsLog.open (o.regionsLog);
    regionsLog << "@This file list merge regions ', ' separated and reads newline separated. \n@The next list separated by an additonal newline\n";
    }


    for(int i = 0; i < recordtable.size(); ++i){
        if(logregions)
            regionsLog << toCString(std::get<3>(table[i])) << "\n";

        for(int j = 0; j < recordtable[i].size(); ++j){
            auto & record = recordtable[i][j];
            if(length(record.seq) == 0){
                ++empty_dups;
//                 std::cout << "emptyDups: " << record.qName << "\n";
            }
            else
            {
                IDMAP::iterator it;
                if (o.rmMulti){
                    it = idMap_occ.find(record.qName);
                    if(it != idMap_occ.end()){
                        if(it->second == 1){
                            uniqueExtractedReads.push_back(std::move(recordtable[i][j]));
                            if(logregions){
                                regionsLog << toCString(recordtable[i][j].qName) << "\n";
                            }
                        }else{
//                             std::cout << "Dups: " << record.qName << "\n";
                            ++dups;
                        }
                    }
                    else
                    {
                        std::cout << "This should not happen\n";
                    }
                }
                else
                {
                    it = idMap.find(record.qName);
                    if(it == idMap.end())
                    {
                        idMap[record.qName] = 1;
                        uniqueExtractedReads.push_back(std::move(recordtable[i][j]));
                        if(logregions){
                            regionsLog << toCString(recordtable[i][j].qName) << "\n";
                        }
                    }
                    else
                    {
//                         std::cout << "Dups: " << record.qName << "\n";
                        ++dups;
                    }
                }
            }
        }
        if(logregions)
            regionsLog << "\n";
    }
    if(logregions)
        regionsLog.close();

    std::cout << "Finished extraction. Found " << uniqueExtractedReads.size() << " unique reads. Removed " << dups << " duplicates. Reads with no sequence: " << empty_dups << "(due to being multimappers)\n";

    return uniqueExtractedReads;
}
