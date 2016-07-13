/* align_seq2_raw.cpp
 *
 * Copyright (C) 2009-2014 Nathan Clement
 *
 * This file is part of GNUMAP.
 *
 * GNUMAP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GNUMAP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNUMAP.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

bool process_hits( Genome& gen, seq_map& unique, const Read& search,
        map< pair< unsigned long, unsigned int >, int >& possible_locs, double& denominator,
        double min_align_score, double& top_align_score, bin_seq& bs,
        int strand, int thread_id )
{
    map< pair< unsigned long, unsigned int >, int >::iterator loc_it;
    for( loc_it = possible_locs.begin(); loc_it != possible_locs.end(); ++loc_it )
    {
        // Define the number of hits that have to be matching at each genomic position
        // before we'll even match it
        if( loc_it->second < gMIN_JUMP_MATCHES )
        {
            continue;
        }

        if( loc_it->second == -1 )
        {
            // we've already aligned this position
            continue;
        }

        unsigned long sa_val = ( gen.get_sa_coord( loc_it->first.first ) <= loc_it->first.second ) ? 0 : ( gen.get_sa_coord( loc_it->first.first ) - loc_it->first.second );
        
        string to_match = gen.GetString( sa_val, search.length );
        
        if( to_match.size() == 0 )
        {	//means it's on a chromosome boundary--not valid sequence
#ifdef DEBUG
            fprintf( stderr, "On chromosome boundary...skipping\n" );
#endif
            continue;
        }

        double align_score = 0;
        char aligned[ to_match.length() + 1 ];
        strcpy( aligned, to_match.c_str() );

        //if(gBISULFITE)
            //align_score = bs.get_align_score_w_traceback(search,to_match,aligned);
        //else
        
        
        if( gNW )
        {
            align_score = bs.get_align_score( search, to_match );
#ifdef DEBUG_NW
            // Record the number of nw's this thread performs
            num_nw[ thread_id ]++;
#endif
        }
        else
        {
            // if we're not doing the NW alignment, set the alignment score to
            // the hit count
            align_score = possible_locs[ loc_it->first ];
        }

        // mark the position as visited
        // this should only make a difference when doing NW alignments
        possible_locs[ loc_it->first ] = -1;

#ifdef DEBUG
        if(gVERBOSE > 0)
        {
            cerr << "(" << sa_val << ",";
            cerr << loc_it->second << ") ";
            cerr << "matching: " << to_match ;
            cerr << "  with  : " << consensus;
            cerr << "\tat pos   " << sa_val;
            cerr << "\tand step " << i;
            cerr << "\tand alignment " << align_score;
            cerr << "\tand min " << min_align_score << endl;
        }
#endif
                            
        if( align_score > top_align_score )
        {
            top_align_score = align_score;
        }
            
        string unique_string = to_match;
        
        if( align_score >= min_align_score )
        {
            //fprintf(stderr,"\t**Aligned!\n");

            // take the rev_c before so we don't need to take it inside the ScoredSeq
            // constructor as well.

            ScoredSeq* temp;
            
            if( gSNP )
            {
                temp = new SNPScoredSeq( to_match, align_score, sa_val, strand );
            }
            else if( gBISULFITE || gATOG )
            {
                temp = new BSScoredSeq(to_match,align_score, sa_val, strand );
                //fprintf(stderr,"New BSScoredSeq\n");
            }
            else
            {
                temp = new NormalScoredSeq( to_match, align_score, sa_val, strand );
            }

            if( strand == NEG_STRAND )
            {
                to_match = reverse_comp( to_match );
            }
            
            seq_map::iterator it = unique.find( to_match );
            
            //pair<seq_map::iterator,bool> it_bool = unique.insert(pair<string,ScoredSeq*>(to_match,temp));

            if( it == unique.end() )
            {	// It wasn't found in the set
                unique.insert(pair< string, ScoredSeq* >( to_match, temp ) );
                denominator += exp( align_score );
#ifdef DEBUG
                //fprintf(stderr,"[%10s%c] just added %s at jump %u with size %lu\n",search.name,strand==POS_STRAND ? '+' : '-',gen.GetPos(loc_it->first,strand).c_str(),i,hashes->size);
                //fprintf(stderr,"[%10s%c]     size of unique is now %lu\n",search.name,strand==POS_STRAND?'+':'-',unique.size());
#endif
            }
            else
            { //if(it != unique.end()) 	// It was found in the set
                
                delete temp;

                if( gUNIQUE )
                {	            // Need to return early because this matches to
                                // more than one location
#ifdef DEBUG
                    cerr << "not a unique mapping! returning no mapping!" << endl;
#endif
                    return false;
                }


                if( ( *it ).second->add_spot( sa_val, strand ) )
                {
                    //fprintf(stderr,"[%s] %lu:%d Not already here\n",search.name, sa_val, strand);
                    denominator += exp( align_score );
                }
                //else
                    //fprintf(stderr,"[%s] %lu:%d Already in here, sorry\n",search.name, sa_val, strand );
#ifdef DEBUG
                fprintf( stderr,"[%10s%c]     Didn't add.  size of unique is now %lu\n", search.name, strand == NEG_STRAND ? '+' : '-', unique.size() );
#endif
            }
            
        }
//			else
//				fprintf(stderr,"\tCould not align with align score %f and min_score %f\n",align_score,min_align_score);
    }
    return true;
}

bool align_sequence(Genome &gen, seq_map &unique, const Read &search, const string &consensus, 
				double min_align_score, double &denominator, double &top_align_score,
				int strand, int thread_id) {

	bin_seq bs;
	
	unsigned int i, j, last_kmer_pos = search.length - gMER_SIZE;
    int sa_num_hits = 0;
	
	// A map of genome positions and the number of times they are referenced
	map< pair< unsigned long, unsigned int >, int> possible_locs;

	for( i = 0; i < last_kmer_pos; i += gJUMP_SIZE)
    {
        uint64_t start, end;
        string consensus_piece;

		// If we've run accross a kmer location that is clear (possibly because we've
		// cleared it in a previous step of the algorithm), we want to get another kmer
		// that's right next to it instead of jumping over everything.
		for( j = 0; j + i < last_kmer_pos; j++ )
        {
			consensus_piece = consensus.substr( i + j, gMER_SIZE );
#ifdef DEBUG
			fprintf( stderr, "Consensus piece: >%s<\n", consensus_piece.c_str() );
#endif
            gen.get_sa_int( consensus_piece, &start, &end );
            if( end == 0 && start == 0 )
            {
                //cerr << i + j << " bad " << endl;
                continue;
            }
            else if( gMAX_KMER_SIZE > 0 && end - start + 1 > gMAX_KMER_SIZE )
            {
                //cerr << i + j << " bad2 " << endl;
                continue;
            }
            else
            {
                sa_num_hits = end - start + 1;
                //cerr << i + j << " done searching!\t" << sa_num_hits << endl;
                //cerr << "At location " << j + i << ", " << sa_num_hits << " hits found." << endl;
                break;
            }

#ifdef DEBUG
			cerr << "At location " << j + i << ", no valid locations found." << endl;
#endif
		} // end of for-loop identifying valid index location

		// Increment i by j so we don't look here again.
		i += j;
		
        if( end == 0 && start == 0 )
        {
            break;
        }
        
        if( gMAX_KMER_SIZE > 0 && end - start + 1 > gMAX_KMER_SIZE )
        {
            break;
        }

#ifdef DEBUG
        if( end != 0 && start != 0 )
        {
            cerr << "At location " << i << ", found " << ( end - start + 1 ) << " hits for this string." << endl;
        }
        else
        {
            cerr << "At location " << i << ", found " << 0 << " hits for this string." << endl;
        }
#endif
        
        if( gMAX_KMER_SIZE > 0 && end - start + 1 > gMAX_KMER_SIZE )
        {
            cerr << "ERROR: num of kmer hits: " << end - start + 1 << endl;
        }

        // Update genome position counts
        for( unsigned int vit = start; vit <= end; vit++ )
        {
            //Match the sequence to a sequence that has two characters before and after.
            //string to_match = gen.GetString((*vit)-i-2, search.size()+4);
            
            //unsigned long beginning = ( gen.get_sa_coord( vit ) <= i ) ? 0 : ( gen.get_sa_coord( vit ) - i );
            
            // Increment the number of times this occurs
            if( possible_locs[ pair< unsigned long, unsigned int >( vit, i ) ] != -1 )
            {
                possible_locs[ pair< unsigned long, unsigned int >( vit, i ) ]++;
            }
        }

        // @masakistan: if we're not doing needleman wunsh alignments, we only have to update the counts
        // and once the total counts have been calculated we can filter for the positions that have enough
        // hits.
        // so, we continue on if we're not doing NW alignments and then do an end calculation at the bottom
        if( !gNW )
        {
            continue;
        }

        // @masakistan: if we are going to do the needleman-wunsch alignments, we will check as we go
        // for genomic positions that have enough hits to permit a hit to be made.
 
        /*bool process_hits( Genome& gen, seq_map& unique, const Read& search,
            map< unsigned long, int >& possible_locs, double& denominator,
            double min_align_score, double& top_align_score, bin_seq& bs,
            int strand, int thread_id )*/
        bool goon = process_hits( gen, unique, search, possible_locs, denominator, min_align_score, top_align_score, bs, strand, thread_id );

        if( !goon )
        {
            return false;
        }
       
		
		if( unique.size() > gMAX_MATCHES )
        {

//#ifdef DEBUG
            //cerr << "too many matches! returning no mapping! " << search.name << "\t" << unique.size() << endl;
//#endif
			return false;
        }
		
		
		// If we want to do the fast mode, only look at one hit location
		if( gFAST)
        {
			break;
        }
	
	}  // end of for over all the hit positions

    // process all the hits now if we're not doing the NW alignments
    if( !gNW )
    {
        bool goon = process_hits( gen, unique, search, possible_locs, denominator, min_align_score, top_align_score, bs, strand, thread_id );
        
        if( !goon )
        {
            cerr << "Turning off Needleman-Wunsch aligments and doing unique alignments are not possible" << endl;
        }
    }
	
	return true;
}
