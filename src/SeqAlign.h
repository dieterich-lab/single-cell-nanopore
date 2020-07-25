// SeqAlign.h

#ifndef FLEXBAR_SEQALIGN_H
#define FLEXBAR_SEQALIGN_H

#include <assert.h>

tbb::mutex ouputMutex;

template <typename TSeqStr, typename TString, class TAlgorithm>
class SeqAlign {

private:

	typedef AlignResults<TSeqStr> TAlignResults;

	const flexbar::LogAlign    m_log;
	const flexbar::FileFormat  m_format;
	const flexbar::PairOverlap m_poMode;

	const bool m_htrimMaxFirstOnly;

	const std::string m_htrimLeft, m_htrimRight;

	const unsigned int m_htrimMinLength, m_htrimMinLength2, m_htrimMaxLength;

	const float m_htrimErrorRate;

    const bool m_barcodeAlignment;

	const bool m_isBarcoding, m_writeTag, m_umiTags, m_strictRegion, m_addBarcodeAdapter, m_printAlignment, m_logEverything;
	const int m_minLength, m_minOverlap, m_tailLength_extension;
	const float m_errorRate;
	const unsigned int m_bundleSize;
    const int m_barcode_umi_length;
    const int m_keepbp;

	tbb::atomic<unsigned long> m_nPreShortReads, m_modified;
	tbb::concurrent_vector<flexbar::TBar> *m_queries;
	tbb::concurrent_vector<unsigned long> m_rmOverlaps;

	std::ostream *m_out, *umi_out;
    AlignmentResults & res;
/*
    typedef std::map<seqan::CharString, short>  PScore;
    PScore & m_leftTailScores;
    PScore & m_rightTailScores;

    typedef std::tuple<seqan::CharString, seqan::Dna5String, seqan::CharString> TRecord;

    std::vector<TRecord> & m_leftTail;
    std::vector<TRecord> & m_rightTail;*/

	TAlgorithm m_algo;

public:

	SeqAlign(tbb::concurrent_vector<flexbar::TBar> *queries, const Options &o, int minOverlap, float errorRate, const int tailLength, const int match, const int mismatch, const int gapCost, const bool isBarcoding, AlignmentResults & res):

			m_minOverlap(minOverlap),
			m_errorRate(errorRate),
			m_tailLength_extension(tailLength),
			m_barcodeAlignment(o.barcodeAlignment),
			m_isBarcoding(isBarcoding),
			m_umiTags(o.umiTags),
			m_minLength(o.min_readLen),
			m_poMode(o.poMode),
			m_log(o.logAlign),
			m_printAlignment(o.printAlignment),
			m_logEverything(o.logEverything),
			m_format(o.format),
			m_writeTag(o.useRemovalTag),
			m_addBarcodeAdapter(o.addBarcodeAdapter),
			m_strictRegion(! o.relaxRegion),
			m_bundleSize(o.bundleSize),
			m_out(o.out),
			umi_out(o.umi_out),
			m_nPreShortReads(0),
			m_modified(0),
			m_barcode_umi_length(o.barcodeUmiLength),
			m_keepbp(o.keepbpOfAdapter),
            m_htrimRight(o.htrimRight),
            m_htrimMinLength(o.htrimMinLength),
            m_htrimMinLength2(o.htrimMinLength2),
            m_htrimMaxLength(o.htrimMaxLength),
            m_htrimMaxFirstOnly(o.htrimMaxFirstOnly),
            m_htrimErrorRate(o.h_errorRate),
            res(res),
			m_algo(TAlgorithm(o, match, mismatch, gapCost, ! isBarcoding)){

		m_queries    = queries;
		m_rmOverlaps = tbb::concurrent_vector<unsigned long>(flexbar::MAX_READLENGTH + 1, 0);
	};
	
	int alignSeqRead(flexbar::TSeqRead* sr, const bool performRemoval, flexbar::Alignments &alignments, flexbar::ComputeCycle &cycle, unsigned int &idxAl, const flexbar::AlignmentMode &alMode, const flexbar::TrimEnd trimEnd, const TSeqStr &addBarcode){

		using namespace std;
		using namespace flexbar;

		using seqan::prefix;
		using seqan::suffix;

		TSeqRead &seqRead = *sr;

		int readLength    = length(seqRead.seq);

		if(! m_isBarcoding && readLength < m_minLength){
			if(cycle != PRELOAD) ++m_nPreShortReads;
			// return 0;
		}

		if(readLength < 1) return 0;

		if(cycle == PRELOAD){

			if(idxAl == 0) reserve(alignments.aset, m_bundleSize * m_queries->size());

			for(unsigned int i = 0; i < m_queries->size(); ++i){

				if     (alMode == ALIGNRCOFF &&   m_queries->at(i).rcAdapter) continue;
				else if(alMode == ALIGNRC    && ! m_queries->at(i).rcAdapter) continue;

				TSeqStr *qseq = &m_queries->at(i).seq;
				TSeqStr *rseq = &seqRead.seq;
				TSeqStr tmp, tmpq;

//                 std::cout << seqRead.id << "\t" << seqRead.seq << "\n";

				if(! m_isBarcoding && m_addBarcodeAdapter && addBarcode != ""){
					tmpq = addBarcode;
					append(tmpq, m_queries->at(i).seq);
					qseq = &tmpq;
				}

				if(trimEnd == LTAIL || trimEnd == RTAIL){
					int tailLength  = length(*qseq) + m_keepbp + m_tailLength_extension;

					if(tailLength < readLength){
						if(trimEnd == LTAIL) tmp = prefix(seqRead.seq, tailLength);
						else                 tmp = suffix(seqRead.seq, readLength - tailLength);
						rseq = &tmp;
					}
				}

				TAlign align;
				appendValue(alignments.aset, align);
				resize(rows(alignments.aset[idxAl]), 2);

				assignSource(row(alignments.aset[idxAl], 0), *rseq);
				assignSource(row(alignments.aset[idxAl], 1), *qseq);

				++idxAl;
			}
			return 0;
		}

		TAlignResults am;

		int qIndex  = -1;
		int pos_bestScore = -1;
		int amScore = numeric_limits<int>::min();
		int amScore_bestScore = amScore;

		std::vector<TAlignResults> am_v;
		std::vector<int> qIndex_v;
		std::vector<int> scores;

		// align each query sequence and store best one
		for(unsigned int i = 0; i < m_queries->size(); ++i){

			if     (alMode == ALIGNRCOFF &&   m_queries->at(i).rcAdapter) continue;
			else if(alMode == ALIGNRC    && ! m_queries->at(i).rcAdapter) continue;

			TAlignResults a;

			// global sequence alignment
			m_algo.alignGlobal(a, alignments, cycle, idxAl++, trimEnd);

			a.queryLength = length(m_queries->at(i).seq);

			if(! m_isBarcoding && m_addBarcodeAdapter && addBarcode != ""){
				a.queryLength += length(addBarcode);
			}

			a.tailLength  = m_tailLength_extension + a.queryLength + m_keepbp;

			a.overlapLength = a.endPos - a.startPos;
			a.allowedErrors = m_errorRate * a.overlapLength;

			float madeErrors = static_cast<float>(a.mismatches + a.gapsR + a.gapsA);
			int minOverlap   = (m_isBarcoding && m_minOverlap == 0) ? a.queryLength : m_minOverlap;

			if(! m_isBarcoding && m_poMode == PON && seqRead.pairOverlap &&
				(trimEnd == RIGHT || trimEnd == RTAIL)) minOverlap = 1;

			bool validAl = true;

			if(((trimEnd == RTAIL  || trimEnd == RIGHT) && a.startPosA < a.startPosS && m_strictRegion) ||
			   ((trimEnd == LTAIL  || trimEnd == LEFT)  && a.endPosA   > a.endPosS   && m_strictRegion) ||
			     a.overlapLength < 1){

				validAl = false;
			}

			// check if alignment is valid, score max, number of errors and overlap length
			if(validAl && a.score > amScore && madeErrors <= a.allowedErrors && a.overlapLength >= minOverlap){
				amScore = a.score;
				if(m_logEverything){
					scores.push_back(amScore);
					am_v.push_back(a);
					qIndex_v.push_back(i);
				}

				if(amScore_bestScore < amScore){
					if(!m_logEverything){
						am      = a;
						qIndex  = i;
					}
					amScore_bestScore = amScore;
					pos_bestScore = qIndex_v.size() - 1;
				}
			}
		}

		// If we are only interested in the best alignment put them in the vector after checking all queries
		//TODO add this to other version
		if(!m_logEverything && qIndex != -1){
			am_v.push_back(am);
			qIndex_v.push_back(qIndex);
		}
		else if (qIndex_v.size() > 1 && am_v.size() > 1)
		{
//			qIndex_v.push_back(qIndex_v[pos_bestScore]);
//			qIndex_v.erase(qIndex_v.begin() + pos_bestScore);
//			am_v.push_back(am_v[pos_bestScore]);
//			am_v.erase(am_v.begin() + pos_bestScore);
			std::iter_swap(qIndex_v.begin() + pos_bestScore, qIndex_v.begin() + (qIndex_v.size() - 1)); // + (qIndex_v.size() - 1)
			std::iter_swap(am_v.begin() + pos_bestScore, am_v.begin() + (am_v.size() - 1));
		}

		int smallest_diff = -1;
		if(m_logEverything){
			sort(scores.rbegin(), scores.rend());
		if(scores.size() > 1)
			smallest_diff = amScore_bestScore - scores[1];
		}

		string smallest_diff_to_best_score;
		if(smallest_diff == -1){
			smallest_diff_to_best_score = "";
		}else{
			smallest_diff_to_best_score = "/" + std::to_string(smallest_diff);
		}

		stringstream s;

        int polyTlength = 0;
        int prefixPolyT = 0;
        bool valid_read = false;

        if(m_htrimRight != ""){

            auto tmpseq = seqRead.seq;
            int nonPolyTprefixLength = m_barcode_umi_length + m_keepbp - 1;
            if(trimEnd == TrimEnd::RTAIL || trimEnd == TrimEnd::RIGHT)
                seqan::reverseComplement(tmpseq);
            valid_read = false;

            for(unsigned int pos = 0; pos < m_htrimRight.length(); ++pos){

                char nuc = m_htrimRight[pos];

//                    unsigned int seqLen = length(seqRead->seq);
                int cutPos = 0;
                int notNuc = 0;

                for(int i = nonPolyTprefixLength; i < readLength; ++i){
                    if(tmpseq[i] != nuc){
                        notNuc++;
                    }
                    else if((notNuc) <= m_htrimErrorRate * (i - nonPolyTprefixLength + 1)){

                        if(m_htrimMaxLength != 0 && (i - nonPolyTprefixLength + 1) > m_htrimMaxLength && (!m_htrimMaxFirstOnly || pos == 0))
                            break;

                        cutPos = i;
                    }

                }

//                 s << "\ncutoff: " << (readLength - cutPos + nonPolyTprefixLength) << "\n";
                unsigned int htrimMinLength = m_htrimMinLength;
                if(m_htrimMinLength2 > 0 && pos > 0) htrimMinLength = m_htrimMinLength2;


//                 s << readLength << "\t" << htrimMinLength << "\t" << cutPos << "\t" << nonPolyTprefixLength << "\n";
                if(cutPos > 0 && (cutPos - nonPolyTprefixLength + 1) >= htrimMinLength){
//                         erase(seqRead->seq, cutPos, length(seqRead->seq));
                    valid_read = true;
                    polyTlength = cutPos;

                    //TODO add option
                    for(int i = nonPolyTprefixLength - 1; i > 0; --i){
                        if(tmpseq[i] == nuc){
                            polyTlength++;
                            prefixPolyT++;
                        }
                        else
                        {
                            break;
                        }
                    }
                }
            }
        }

        //TODO m_queries->at(qIndex).rcAdapter does not work but from qIndex_v[i] rc can be derived
        typedef std::map<seqan::CharString, uint32_t>  PScore;
		// valid alignment
		if(qIndex_v.size() > 0){
			for(int i = 0; i < qIndex_v.size(); ++i){
//                 if(!valid_read)
//                     s << "Not valid; No oligoT found\n";
				TSeqRead seqReadTmp = seqRead;
                TSeqRead barcode = seqRead;
				if(!m_logEverything)
				{
					//only do one iteration of the loop
	//				i = qIndex_v.size();
					// use best alignment
					qIndex = qIndex_v.front();
					am = am_v.front();
				}
				else
				{
					if(i < qIndex_v.size() - 1){
// 						if(m_log == ALL)
// 							s << "Alternative alignment:\n";
						qIndex = qIndex_v[i];
						am = am_v[i];
					}
					else
					{
// 						if(m_log == ALL)
// 							s << "Best alignment (" << qIndex_v.size() << smallest_diff_to_best_score << "):" << "\n";
						qIndex = qIndex_v[i];
						am = am_v[i];
					}
				}
				TrimEnd trEnd = trimEnd;

				// trim read based on alignment
				if(performRemoval){

					if(trEnd == ANY){
						if(am.startPosA <= am.startPosS && am.endPosS <= am.endPosA){
							seqReadTmp.seq = "";
							if(m_format == FASTQ) seqReadTmp.qual = "";
						}
						else if(qIndex_v[i] == 1){
							trEnd = RIGHT;
						}
						else trEnd = LEFT;
					}

					switch(trEnd){

						int rCutPos;

						case LTAIL:
						case LEFT:
                        {
							rCutPos = am.endPos;

							// translate alignment end pos to read idx
							if(am.startPosS > 0) rCutPos -= am.startPosS;

							// adjust to inner read gaps
							rCutPos -= am.gapsR;

							if(rCutPos > readLength) rCutPos = readLength;

                            int modCut = (rCutPos > m_keepbp) ? (rCutPos - m_keepbp) : 0;
                            if(!m_barcodeAlignment){
                                erase(seqReadTmp.seq, 0, modCut);
                                erase(barcode.seq, rCutPos, length(barcode.seq));
                            }
                            else
                            {
                                erase(seqReadTmp.seq, 0, modCut);
                                erase(barcode.seq, rCutPos, length(barcode.seq));
                            }


							if(m_format == FASTQ){
								erase(seqReadTmp.qual, 0, rCutPos);
                                erase(barcode.qual, rCutPos, length(barcode.qual));
                            }
                            //TODO add modified removal to fastq

							break;
                        }

						case RTAIL:
							// adjust cut pos to original read length
							am.startPos += readLength - am.tailLength;

						case RIGHT:
                        {
							rCutPos = am.startPos;

							// skipped restriction
							if(rCutPos < 0) rCutPos = 0;

                            int modCut2 = (rCutPos + m_keepbp < readLength) ? (rCutPos + m_keepbp) : (readLength - 1);
                            if(!m_barcodeAlignment){
                                erase(seqReadTmp.seq, modCut2, readLength);
                                erase(barcode.seq, 0, modCut2);
                            }
                            else
                            {
                                erase(seqReadTmp.seq, rCutPos, readLength);
                                erase(barcode.seq, 0, rCutPos);
                            }

							if(m_format == FASTQ){
								erase(seqReadTmp.qual, rCutPos, readLength);
                                erase(barcode.qual, 0, rCutPos);
                                assert(length(seqReadTmp.seq) == length(seqReadTmp.qual));
                            }

							break;
                        }

			     			case ANY:;
					}

					++m_modified;
					if(! m_isBarcoding){
						if(! m_queries->at(qIndex).rcAdapter){
							seqReadTmp.rmAdapter   = true;
						}else{
		                                  seqReadTmp.rmAdapterRC = true;
						}
					}
					// count number of removals for each query
					m_queries->at(qIndex).rmOverlap++;
					if(am.overlapLength == am.queryLength)
					m_queries->at(qIndex).rmFull++;
					if(m_writeTag){
						append(seqReadTmp.id, "_Flexbar_removal");

						if(! m_isBarcoding){
							append(seqReadTmp.id, "_");
							append(seqReadTmp.id, m_queries->at(qIndex).id);
						}
					}
					// store overlap occurrences
					if(am.overlapLength <= MAX_READLENGTH) m_rmOverlaps.at(am.overlapLength)++;
					else cerr << "\nCompile Flexbar with larger max read length for correct overlap stats.\n" << endl;
				}

				// valid alignment, not neccesarily removal
				if(m_umiTags && am.umiTag != ""){
					append(seqReadTmp.umi, "_");
					append(seqReadTmp.umi, am.umiTag);
				}

				// alignment stats
				if(m_log == ALL || (m_log == MOD && performRemoval)){
 					if(performRemoval){
 						s << "Sequence removal:";

 						     if(trEnd == LEFT  || trEnd == LTAIL) s << " left side\n";
 						else if(trEnd == RIGHT || trEnd == RTAIL) s << " right side\n";
 						else                                      s << " any side\n";
 					}
 					else s << "Sequence detection, no removal:\n";

					s << "  query id         " << m_queries->at(qIndex).id            << "\n"
	  				  << "  query pos        " << am.startPosA << "-" << am.endPosA   << "\n"
					  << "  read id          " << seqReadTmp.id                          << "\n"
					  << "  read pos         " << am.startPosS << "-" << am.endPosS   << "\n"
					  << "  score            " << am.score                            << "\n"
					  << "  overlap          " << am.overlapLength                    << "\n"
					  << "  errors           " << am.gapsR + am.gapsA + am.mismatches << "\n"
					  << "  error threshold  " << am.allowedErrors                    << "\n";

				if(performRemoval){
						s << "  remaining read   " << seqReadTmp.seq << "\n";

					if(m_format == FASTQ)
						s << "  remaining qual   " << seqReadTmp.qual << "\n";
				}
                if(performRemoval){
						s << "  barcode seq      " << barcode.seq << "\n";

					if(m_format == FASTQ)
						s << "  barcode qual     " << barcode.qual << "\n";
				}
				s << "\n  Alignment:\n" << endl << am.alString;


                }
			else if(m_log == TAB && m_barcodeAlignment){

                if(m_printAlignment)
                    s << am.alString;
                s << seqReadTmp.id   << "\t";

                seqan::Dna5String barcode = m_queries->at(qIndex).seq;
                if(trimEnd == RTAIL  || trimEnd == RIGHT)
                {
                    seqan::reverseComplement(barcode);
                    PScore::iterator it = res.rightTailScores.find(seqReadTmp.id);
                    if(it != res.rightTailScores.end())
                    {

                        s << m_queries->at(qIndex).id << "\tno\tright\t" << it->second << "\t";
                    }
                    else
                    {
                        s << m_queries->at(qIndex).id << "\tno\tright\tWARNING VALUE NOT FOUND\t";
                        std::cerr << "WARNING NO PRIMER ALIGNMENT SCORE FOR READID" << seqReadTmp.id << "FOUND\n";
                    }
                }

                if(trimEnd == LTAIL  || trimEnd == LEFT){
                    PScore::iterator it = res.leftTailScores.find(seqReadTmp.id);
                    if(it != res.leftTailScores.end())
                    {
                        s << m_queries->at(qIndex).id << "\tno\tleft\t" << it->second << "\t";
                    }
                    else
                    {
                        s << m_queries->at(qIndex).id << "\tno\tleft\tWARNING VALUE NOT FOUND\t";
                        std::cerr << "WARNING NO PRIMER ALIGNMENT SCORE FOR READID" << seqReadTmp.id << "FOUND\n";
                    }
                }

                //adaptor_score
                int barcodeStartOffset = (trimEnd == RTAIL  || trimEnd == RIGHT) ? (am.endPosS - am.endPosA) : am.startPosA;
                int barcodeLength = am.endPosA - am.startPosA + am.gapsR;
		int umiLength = m_barcode_umi_length - prefixPolyT - barcodeLength;

//                 std::cout << am.gapsR  << "\t" << am.gapsA << "\n";
//                 std::cout << barcodeLength << "\t" << am.startPosS  << "\t" << am.endPosS  << "\t" << am.startPosA << "\t" << am.endPosA << "\n";

                s << am.score << "\t" << am.gapsR + am.gapsA  << "\t" << (barcodeStartOffset - m_keepbp) << "\t" << am.mismatches << "\t";
                if(valid_read)
                    s << umiLength << "\t" << polyTlength - m_barcode_umi_length << "\t";
                else
                    s << length(seqReadTmp.seq) << "\t" << "-1" << "\t";
                if(i < qIndex_v.size() - 1)
                    s << "yes";
                else
                    s << "no";
                
                //if(length(seqReadTmp.seq) >= m_keepbp + 15)
                umi_out << ">" << seqReadTmp.id << "\t" << barcode << std::endl << seqReadTmp.seq << std::endl;
                
		s << "\t" << "0\n";
            }

				if(i == qIndex_v.size() - 1 || !m_logEverything){
					seqRead = seqReadTmp;
				}


				//TODO new mutexLock!!!!
// 				std::cout << "Primer alignment: " << m_barcodeAlignment << "\n";


                if(!m_barcodeAlignment){
                    ouputMutex.lock();
                    uint32_t alignmentScore = (am.score > 0) ? am.score : 0;
                    if(length(seqReadTmp.seq) <= (m_barcode_umi_length - 10))
                    {
                        s << seqRead.id << "\t" << "NA\t" << "no\t" << alignmentScore << "\t-1\t" << "-1\t" <<  "-1\t" <<  "-1\t" <<  "-1\t" << "-1\t" << "-1\t" <<  "no\t" << "0\n";
                    }
                    else if(trEnd == LEFT  || trEnd == LTAIL)
                    {
                        res.leftTail.push_back(make_tuple(seqReadTmp.id, seqReadTmp.seq, seqReadTmp.qual));
                        res.leftTailScores[static_cast<seqan::CharString>(seqReadTmp.id)] = alignmentScore;
                    }
                    else if(trEnd == RIGHT || trEnd == RTAIL) //TODO think about this value
                    {
                        res.rightTail.push_back(make_tuple(seqReadTmp.id, seqReadTmp.seq, seqReadTmp.qual));
                        res.rightTailScores[static_cast<seqan::CharString>(seqReadTmp.id)] = alignmentScore;
                    }
                    ouputMutex.unlock();
                }

			}
		}
		else
        {
            //TODO add unvalid
//             unvalid.push_back(make_tuple(seqRead.id, seqRead.seq, seqRead.qual));
            if(m_log == ALL){
                s << "Unvalid alignment:"        << "\n"
                << "read id   " << seqRead.id  << "\n"
                << "read seq  " << seqRead.seq << "\n\n" << endl;
            }

            if(m_log == TAB){
//                 <read_id>,<barcode>,<adaptor_invalid>,<adaptor_side>,<adaptor_score>,<barcode_score>,<barcode_indel>,<barcode_start>,<barcode_mismatch>,<umi_length>,<polyT_length>,<barcode_alt>
                if(!m_barcodeAlignment){
                    s << seqRead.id << "\t" << "NA\t" << "yes\t" << "-1\t" << "-1\t" << "-1\t" <<  "-1\t" <<  "-1\t" <<  "-1\t" << "-1\t" << "-1\t" <<  "no\t" << "0\n";
                }
                else
                {
                    s << seqRead.id << "\t" << "NA\t" << "no\t";

                    if(trimEnd == RTAIL  || trimEnd == RIGHT)
                    {
                        PScore::iterator it = res.rightTailScores.find(seqRead.id);
                        if(it != res.rightTailScores.end())
                        {
                            s << "right\t" << it->second << "\t";
                        }
                        else
                        {
                            s << "right\t" << "WARNING VALUE NOT FOUND" << "\t";
                            std::cerr << "WARNING NO PRIMER ALIGNMENT SCORE FOR READID" << seqRead.id << "FOUND\n";
                        }
                    }

                    if(trimEnd == LTAIL  || trimEnd == LEFT){
                        PScore::iterator it = res.leftTailScores.find(seqRead.id);
                        if(it != res.leftTailScores.end())
                        {
                            s << "left\t" << it->second << "\t";
                        }
                        else
                        {
                            s << "left\t" << "WARNING VALUE NOT FOUND" << "\t";
                            std::cerr << "WARNING NO PRIMER ALIGNMENT SCORE FOR READID" << seqRead.id << "FOUND\n";
                        }
                    }

                    s << "-1\t" <<  "-1\t" <<  "-1\t" << "-1\t" << "-1\t" << "-1\t" << "-1\t" <<  "no\t" << "0\n";
                }

            }

		}

		ouputMutex.lock();
		*m_out << s.str();
		ouputMutex.unlock();

		return ++qIndex;
	}

	std::string getOverlapStatsString(){

		using namespace std;
		using namespace flexbar;

		unsigned long nValues = 0, halfValues = 0, cumValues = 0, lenSum = 0;
		unsigned int max = 0, median = 0, mean = 0;

		unsigned int min = numeric_limits<unsigned int>::max();

		for(unsigned int i = 0; i <= MAX_READLENGTH; ++i){
			unsigned long lenCount = m_rmOverlaps.at(i);

			if(lenCount > 0 && i < min) min = i;
			if(lenCount > 0 && i > max) max = i;

			nValues += lenCount;
			lenSum  += lenCount * i;
		}

		halfValues = nValues / 2;

		for(unsigned int i = 0; i <= MAX_READLENGTH; ++i){
			cumValues += m_rmOverlaps.at(i);

			if(cumValues >= halfValues){
				median = i;
				break;
			}
		}

		if(m_modified > 0) mean = lenSum / m_modified;

		stringstream s;

		s << "Min, max, mean and median overlap: ";
		s << min << " / " << max << " / " << mean << " / " << median;

		return s.str();
	}

	unsigned long getNrPreShortReads() const {
		return m_nPreShortReads;
	}

	unsigned long getNrModifiedReads() const {
		return m_modified;
	}

};

#endif
