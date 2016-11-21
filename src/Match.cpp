#include "Match.h"
#include "GERMLINE.h"

void Match::extendBack()
{
	// save old position values
	unsigned int SAVE_pms = position_ms;
	
	// iterate backwards through genome
	while(position_ms > 0)
	{
		position_ms--;
		if( !approxEqual() )
		{
			position_ms++;
			break;
		}
	}
	start_ms = position_ms;
	// restore saved values
	position_ms = SAVE_pms;
}

bool Match::approxEqual()
{

		
	// homozygosity check
	if ( node[0] == node[1] )
	{
		if ( ALLOW_HOM )
		{
			if ( (int) ( node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 1 )->getMarkerSet()->getMarkerBits() ).count() 
				 <= ( MAX_ERR_HOM + MAX_ERR_HET ) ) return true; else return false;
		}
		else
		{
			return false;
		}
	}
	else
	{
		// 1. Haplotype extension

		 if ( HAPLOID )
                {
                                if ( (int)(node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet()->getMarkerBits()).count() <= MAX_ERR_HOM ) return true;

		} else
		{
			for ( int a = 0 ; a < 2 ; a++ ) {
				for ( int b = 0 ; b < 2 ; b++ ) { 
					if ( (int)(node[0]->getChromosome( a )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( b )->getMarkerSet()->getMarkerBits()).count() <= MAX_ERR_HOM )
					{
						return true;
					}
				}
			}
		}

		if ( HAPLOID || HAP_EXT ) return false;

		// 2. Genotype extension
		// identify common homozygous SNPs
		boost::dynamic_bitset<> mask
			= ( node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[0]->getChromosome( 1 )->getMarkerSet()->getMarkerBits() ).flip()
			& ( node[1]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 1 )->getMarkerSet()->getMarkerBits() ).flip();

		// assert that homozygous SNPs are identical
		if ( (int) ((node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet()->getMarkerBits()) & mask).count() <= MAX_ERR_HET )
		{
			return true;
		}
		else return false;
	}
}

int Match::scanLeft( unsigned int ms )
{
	bool err = false;
	int marker = MARKER_SET_SIZE - 1;

	if ( HAPLOID ) {
		for ( marker = MARKER_SET_SIZE - 1 ; marker >= 0 && ! err; marker-- )
			if ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits()[marker] != node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits()[marker] )
				err = true;
	}

        else if(node[0]==node[1])
        {
                int cur_marker;
                err=false;
        if ( ALLOW_HOM )
        {
                        int count=0;
                       for (  cur_marker = MARKER_SET_SIZE - 1 ; cur_marker >0 && count<=ERR_W; cur_marker-- )
                                {
                                        if ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits()[cur_marker] != node[1]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits()[cur_marker] )
                                                count++;

                                }
                                if ( cur_marker > marker ) marker = cur_marker;

        }
        }


	 else if ( HAP_EXT )
	{
		int cur_marker;
		for ( int a = 0 ; a < 2 ; a++ ) {
			for ( int b = 0 ; b < 2 ; b++ ) { 
				err = false;
				for ( cur_marker = MARKER_SET_SIZE - 1 ; cur_marker >= 0 && ! err; cur_marker-- )
				{
					if ( node[0]->getChromosome( a )->getMarkerSet(ms)->getMarkerBits()[cur_marker] != node[1]->getChromosome( b )->getMarkerSet(ms)->getMarkerBits()[cur_marker] )
						err = true;
				}
				if ( cur_marker < marker ) marker = cur_marker;
			}
		}
	}

	 else
	{
	boost::dynamic_bitset<> mask
		= ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[0]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip()
		& ( node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip();
	mask = ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits()) & mask;

	for( marker = MARKER_SET_SIZE - 1 ; marker >= 0 && !err ; marker-- )
		if ( mask[marker] ) err = true;
	}
	
	return marker;
}

int Match::scanRight( unsigned int ms )
{
	bool err = false;
	int marker = 0;

	if ( HAPLOID ) {
		for ( marker = 0 ; marker < MARKER_SET_SIZE && ! err; marker++ )
			if ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits()[marker] != node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits()[marker] )
				err = true;
	}

	 else if(node[0]==node[1])
        {
                int cur_marker;
                err=false;
        if ( ALLOW_HOM )
        {
                        int count=0;
                       for (  cur_marker = 0 ; cur_marker < MARKER_SET_SIZE && count<=ERR_W; cur_marker++ )
                                {
                                        if ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits()[cur_marker] != node[1]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits()[cur_marker] )
                                                count++;

                                }
                                if ( cur_marker > marker ) marker = cur_marker;

        }
        }
	 else if ( HAP_EXT )
	{
		int cur_marker;
		for ( int a = 0 ; a < 2 ; a++ ) {
			for ( int b = 0 ; b < 2 ; b++ ) { 
				err = false;
				for ( cur_marker = 0 ; cur_marker < MARKER_SET_SIZE && ! err; cur_marker++ )
				{
					if ( node[0]->getChromosome( a )->getMarkerSet(ms)->getMarkerBits()[cur_marker] != node[1]->getChromosome( b )->getMarkerSet(ms)->getMarkerBits()[cur_marker] )
						err = true;
				}
				if ( cur_marker > marker ) marker = cur_marker;
			}
		}
	}
	else
	{
	boost::dynamic_bitset<> mask
		= ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[0]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip()
		& ( node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip();
	mask = ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits()) & mask;

	for( marker = 0 ; marker < MARKER_SET_SIZE && !err ; marker++ )
		if ( mask[marker] ) err = true;
	}	
	return marker;
}

int Match::diff( unsigned int ms )
{
	boost::dynamic_bitset<> mask
		= ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[0]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip()
		& ( node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip();
	mask = ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits()) & mask;

	return int(mask.count());
}

bool Match::isHom( int n , unsigned int ms )
{
	return (int) ( node[n]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[n]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).count() <= ( MAX_ERR_HOM + MAX_ERR_HET );
}

void Match::print( ostream& fout )
{
	// extend this match from both ends
	unsigned int snp_start = ALL_SNPS.getROIStart().getMarkerNumber() + start_ms * MARKER_SET_SIZE;
	unsigned int snp_end = ALL_SNPS.getROIStart().getMarkerNumber() + ( end_ms + 1 ) * MARKER_SET_SIZE - 1;
	int marker;

	
	if ( WIN_EXT )
	{
		// backwards
		if( start_ms > 0 )
		{
			marker = scanLeft( start_ms - 1 );
			snp_start -= (MARKER_SET_SIZE - marker - 2);
		}
	}
	if ( WIN_EXT || end_ms == num_sets - 2 )
	{
		// forwards
		if( end_ms < num_sets - 1 )
		{
			marker = scanRight( end_ms + 1 );
			snp_end += marker - 1;
		}
	}
	

	bool genetic;
	float distance;
	if ( ( distance = ALL_SNPS.getDistance(snp_start,snp_end,genetic)) < MIN_MATCH_LEN ) return;
	// print

	// get hamming distance & ignored bit count
	int dif = 0;
	for( unsigned int i = start_ms; i <= end_ms ; i++) { dif += diff( i ); }
	
	// calculate if homozygous
	bool hom[2];
	if ( node[0] == node[1] ) { hom[0] = hom[1] = 1; }
	else
	{
		for ( int n = 0 ; n < 2 ; n++ )
		{
			hom[n] = true;
			for ( unsigned int i = start_ms ; i<= end_ms && hom ; i++ )
			{
				hom[n] = isHom( n , i );
			}
		}
	}

	if ( BINARY_OUT )
	{
		unsigned int pid[2];
		pid[0] = node[0]->getNumericID();
		pid[1] = node[1]->getNumericID();
		unsigned int sid[2];
		sid[0] = ALL_SNPS.getSNP(snp_start).getMarkerNumber();
		sid[1] = ALL_SNPS.getSNP(snp_end).getMarkerNumber();
		fout.write( (char*) &pid[0] , sizeof( unsigned int ) );
		fout.write( (char*) &pid[1] , sizeof( unsigned int ) );
		fout.write( (char*) &sid[0] , sizeof( unsigned int ) );
		fout.write( (char*) &sid[1] , sizeof( unsigned int ) );
		fout.write( (char*) &dif , sizeof( int ) );
		fout.write( (char*) &hom[0] , sizeof( bool ) );
		fout.write( (char*) &hom[1] , sizeof( bool ) );
		if(HAPLOID&&node[0]!=node[1])
		{

			string sub1=node[0]->getID(), sub2=node[1]->getID();
			int n1=sub1.find(".0 ");if(n1<0)n1=sub1.find(".1 ");
			int n2=sub2.find(".0 ");if(n2<0)n2=sub2.find(".1 ");
			if(n1>0&&n2>0){
			sub1=sub1.substr(0,n1-2);
			sub2=sub2.substr(0,n2-2);
			}
			if(sub1.compare(sub2)==0)
			{
				                unsigned int pid[2];
                pid[0] = node[0]->getNumericID();
                pid[1] = node[1]->getNumericID();
                unsigned int sid[2];
                sid[0] = ALL_SNPS.getSNP(snp_start).getMarkerNumber();
                sid[1] = ALL_SNPS.getSNP(snp_end).getMarkerNumber();
                fout.write( (char*) &pid[0] , sizeof( unsigned int ) );
                fout.write( (char*) &pid[1] , sizeof( unsigned int ) );
                fout.write( (char*) &sid[0] , sizeof( unsigned int ) );
                fout.write( (char*) &sid[1] , sizeof( unsigned int ) );
                fout.write( (char*) &dif , sizeof( int ) );
                fout.write( (char*) &hom[0] , sizeof( bool ) );
                fout.write( (char*) &hom[1] , sizeof( bool ) );


			}

		}
	}
	else
	//if (IBD )
	{
		stringstream ss;
		string eachline="";
	
		if(!REDUCE){	
			MATCH_FILE << node[0]->getID() << '\t';
			MATCH_FILE << node[1]->getID() << '\t';
		}
		else
		{
			string str;
			int n1=node[0]->getID().find(" "), n2=node[0]->getID().length();
//			if(!HAPLOID)
				/*str=node[0]->getID().substr(0,n1-2);*/
			str=node[0]->getID();
//			else
//				str=node[0]->getID().substr(n1,n2-n1-2);
			ss<<str<<'\t';		//MATCH_FILE2<<str<<'\t';
			n1=node[1]->getID().find(" "); n2=node[1]->getID().length();
  //                      if(!HAPLOID)
                                /*str=node[1]->getID().substr(0,n1-2);*/
								str=node[1]->getID();
    //                    else
       //                         str=node[1]->getID().substr(n1,n2-n1-2);
              ss<<str<<'\t';                  //MATCH_FILE2<<str<<'\t';


		}
		if(!REDUCE)
			MATCH_FILE << ALL_SNPS.getSNP(snp_start).getChr() << '\t';

		
		ss<<ALL_SNPS.getSNP(snp_start).getPhysPos() << '\t';     //MATCH_FILE2 << ALL_SNPS.getSNP(snp_start).getPhysPos() << ' ';
		ss<<ALL_SNPS.getSNP(snp_end).getPhysPos() << '\t';	//MATCH_FILE2 << ALL_SNPS.getSNP(snp_end).getPhysPos() << '\t';
		
		if(!REDUCE){
			MATCH_FILE << ALL_SNPS.getSNP(snp_start).getSNPID() << ' ';
			MATCH_FILE << ALL_SNPS.getSNP(snp_end).getSNPID() << '\t';
		}
		
		ss<<( snp_end - snp_start + 1) << '\t';	//MATCH_FILE2 << ( snp_end - snp_start + 1) << '\t';

		ss<<setiosflags(ios::fixed) << setprecision(2) << distance;	//MATCH_FILE2 << setiosflags(ios::fixed) << setprecision(2) << distance;

		if(!REDUCE){
			MATCH_FILE<<'\t';
			if ( genetic ) MATCH_FILE << "cM" << '\t'; else MATCH_FILE << "MB" << '\t';
		
			MATCH_FILE << dif;
			for ( int n = 0 ; n < 2 ; n++ )
				if ( hom[n] ) MATCH_FILE << '\t' << 1; else MATCH_FILE << '\t' << 0;
		}
		//ss<<endl;	//MATCH_FILE2 << endl;
		ss>>eachline>>eachline;
		MATCH_FILE<<eachline<<'\t';	//person1 2nd col

		ss>>eachline>>eachline;
		MATCH_FILE<<eachline;	//person2 2nd col

		getline(ss,eachline);	//get remaining string from stringstream
		MATCH_FILE<<eachline<<endl;
		ss.str("");
		ss.clear();
	}
	num_matches++;
}

