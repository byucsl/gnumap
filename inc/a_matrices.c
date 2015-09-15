/* a_matrices.c
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
 *
 * This must be separate from the previous method becuase part of command-line parsing
 * looks for the gGAP score.  This will also be changed later if there's a PWM file
 * included
 */
void setup_alignment_matrices() {

/*
	//set up alignment matrix
	for(unsigned int i=0; i<256; i++)
		for(unsigned int j=0; j<4; j++)
			gALIGN_SCORES[i][j] = gGAP;
	gALIGN_SCORES[(int)'a'][0] = gMATCH;
	gALIGN_SCORES[(int)'a'][1] = gMISMATCH;
	gALIGN_SCORES[(int)'a'][2] = gMISMATCH;
	gALIGN_SCORES[(int)'a'][3] = gMISMATCH;

	gALIGN_SCORES[(int)'c'][1] = gMATCH;
	gALIGN_SCORES[(int)'c'][0] = gMISMATCH;
	gALIGN_SCORES[(int)'c'][2] = gMISMATCH;
	gALIGN_SCORES[(int)'c'][3] = gMISMATCH;

	gALIGN_SCORES[(int)'g'][2] = gMATCH;
	gALIGN_SCORES[(int)'g'][0] = gMISMATCH;
	gALIGN_SCORES[(int)'g'][1] = gMISMATCH;
	gALIGN_SCORES[(int)'g'][3] = gMISMATCH;

	gALIGN_SCORES[(int)'t'][3] = gMATCH;
	gALIGN_SCORES[(int)'t'][0] = gMISMATCH;
	gALIGN_SCORES[(int)'t'][1] = gMISMATCH;
	gALIGN_SCORES[(int)'t'][2] = gMISMATCH;
*/

	// These are defined in inc/const_define.h
	// Adjust them so they don't result in overflow with longer sequences
	gMATCH *= gADJUST;
	gTRANSITION *= gADJUST;
	gTRANSVERSION *= gADJUST;
	gGAP *= gADJUST;

	// This is now the DEFAULT scoring system
	//set up alignment matrix
	for(unsigned int i=0; i<256; i++)
		for(unsigned int j=0; j<4; j++)
			gALIGN_SCORES[i][j] = gTRANSVERSION;
	gALIGN_SCORES[(int)'a'][0] = gALIGN_SCORES[(int)'A'][0] = gMATCH;
	gALIGN_SCORES[(int)'a'][1] = gALIGN_SCORES[(int)'A'][1] = gTRANSVERSION;
	gALIGN_SCORES[(int)'a'][2] = gALIGN_SCORES[(int)'A'][2] = gTRANSITION;
	gALIGN_SCORES[(int)'a'][3] = gALIGN_SCORES[(int)'A'][3] = gTRANSVERSION;

	gALIGN_SCORES[(int)'c'][0] = gALIGN_SCORES[(int)'C'][0] = gTRANSVERSION;
	gALIGN_SCORES[(int)'c'][1] = gALIGN_SCORES[(int)'C'][1] = gMATCH;
	gALIGN_SCORES[(int)'c'][2] = gALIGN_SCORES[(int)'C'][2] = gTRANSVERSION;
	gALIGN_SCORES[(int)'c'][3] = gALIGN_SCORES[(int)'C'][3] = gTRANSITION;

	gALIGN_SCORES[(int)'g'][0] = gALIGN_SCORES[(int)'G'][0] = gTRANSITION;
	gALIGN_SCORES[(int)'g'][1] = gALIGN_SCORES[(int)'G'][1] = gTRANSVERSION;
	gALIGN_SCORES[(int)'g'][2] = gALIGN_SCORES[(int)'G'][2] = gMATCH;
	gALIGN_SCORES[(int)'g'][3] = gALIGN_SCORES[(int)'G'][3] = gTRANSVERSION;

	gALIGN_SCORES[(int)'t'][0] = gALIGN_SCORES[(int)'T'][0] = gTRANSVERSION;
	gALIGN_SCORES[(int)'t'][1] = gALIGN_SCORES[(int)'T'][1] = gTRANSITION;
	gALIGN_SCORES[(int)'t'][2] = gALIGN_SCORES[(int)'T'][2] = gTRANSVERSION;
	gALIGN_SCORES[(int)'t'][3] = gALIGN_SCORES[(int)'T'][3] = gMATCH;

	
	/*
	 * set up PHMM alignment matrix
	 * These have been set by Evan.  Eventually, we'll let a flag control this.
	 */
	// This should be the DEFAULT
	
	// Set everything ('n' gPHMM_matching to anything not defined) as 0.005
	for(unsigned int i=0; i<256; i++)
		for(unsigned int j=0; j<4; j++)
			gPHMM_ALIGN_SCORES[i][j] = gPHMM_nsyn;

	gPHMM_ALIGN_SCORES[(int)'a'][0] = gPHMM_ALIGN_SCORES[(int)'A'][0] = gPHMM_match;
	gPHMM_ALIGN_SCORES[(int)'a'][1] = gPHMM_ALIGN_SCORES[(int)'A'][1] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'a'][2] = gPHMM_ALIGN_SCORES[(int)'A'][2] = gPHMM_syn;
	gPHMM_ALIGN_SCORES[(int)'a'][3] = gPHMM_ALIGN_SCORES[(int)'A'][3] = gPHMM_nsyn;

	gPHMM_ALIGN_SCORES[(int)'c'][0] = gPHMM_ALIGN_SCORES[(int)'C'][0] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'c'][1] = gPHMM_ALIGN_SCORES[(int)'C'][1] = gPHMM_match;
	gPHMM_ALIGN_SCORES[(int)'c'][2] = gPHMM_ALIGN_SCORES[(int)'C'][2] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'c'][3] = gPHMM_ALIGN_SCORES[(int)'C'][3] = gPHMM_syn;

	gPHMM_ALIGN_SCORES[(int)'g'][0] = gPHMM_ALIGN_SCORES[(int)'G'][0] = gPHMM_syn;
	gPHMM_ALIGN_SCORES[(int)'g'][1] = gPHMM_ALIGN_SCORES[(int)'G'][1] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'g'][2] = gPHMM_ALIGN_SCORES[(int)'G'][2] = gPHMM_match;
	gPHMM_ALIGN_SCORES[(int)'g'][3] = gPHMM_ALIGN_SCORES[(int)'G'][3] = gPHMM_nsyn;

	gPHMM_ALIGN_SCORES[(int)'t'][0] = gPHMM_ALIGN_SCORES[(int)'T'][0] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'t'][1] = gPHMM_ALIGN_SCORES[(int)'T'][1] = gPHMM_syn;
	gPHMM_ALIGN_SCORES[(int)'t'][2] = gPHMM_ALIGN_SCORES[(int)'T'][2] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'t'][3] = gPHMM_ALIGN_SCORES[(int)'T'][3] = gPHMM_match;

	gPHMM_ALIGN_SCORES[(int)'n'][0] = gPHMM_ALIGN_SCORES[(int)'N'][0] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'n'][1] = gPHMM_ALIGN_SCORES[(int)'N'][1] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'n'][2] = gPHMM_ALIGN_SCORES[(int)'N'][2] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'n'][3] = gPHMM_ALIGN_SCORES[(int)'N'][3] = gPHMM_nsyn;
	
	// This is just to run a test (for now)
	/*
	for(unsigned int i=0; i<256; i++)
		for(unsigned int j=0; j<4; j++)
			gPHMM_ALIGN_SCORES[i][j] = 0.002;

	gPHMM_ALIGN_SCORES[(int)'a'][0] = 0.99;
	gPHMM_ALIGN_SCORES[(int)'a'][1] = 0.002;
	gPHMM_ALIGN_SCORES[(int)'a'][2] = 0.006;
	gPHMM_ALIGN_SCORES[(int)'a'][3] = 0.002;

	gPHMM_ALIGN_SCORES[(int)'c'][0] = 0.002;
	gPHMM_ALIGN_SCORES[(int)'c'][1] = 0.90;
	gPHMM_ALIGN_SCORES[(int)'c'][2] = 0.002;
	gPHMM_ALIGN_SCORES[(int)'c'][3] = 0.006;

	gPHMM_ALIGN_SCORES[(int)'g'][0] = 0.006;
	gPHMM_ALIGN_SCORES[(int)'g'][1] = 0.002;
	gPHMM_ALIGN_SCORES[(int)'g'][2] = 0.99;
	gPHMM_ALIGN_SCORES[(int)'g'][3] = 0.002;

	gPHMM_ALIGN_SCORES[(int)'t'][0] = 0.002;
	gPHMM_ALIGN_SCORES[(int)'t'][1] = 0.006;
	gPHMM_ALIGN_SCORES[(int)'t'][2] = 0.002;
	gPHMM_ALIGN_SCORES[(int)'t'][3] = 0.99;

	gPHMM_ALIGN_SCORES[(int)'n'][0] = 0.002;
	gPHMM_ALIGN_SCORES[(int)'n'][1] = 0.002;
	gPHMM_ALIGN_SCORES[(int)'n'][2] = 0.002;
	gPHMM_ALIGN_SCORES[(int)'n'][3] = 0.002;
	*/

	/*
	float gPHMM_syn=0.025;
	float gPHMM_match=0.95;
	float gPHMM_nsyn=0.0125;
	for(unsigned int i=0; i<256; i++)
		for(unsigned int j=0; j<4; j++)
			gPHMM_ALIGN_SCORES[i][j] = gPHMM_nsyn;

	gPHMM_ALIGN_SCORES[(int)'a'][0] = gPHMM_match;
	gPHMM_ALIGN_SCORES[(int)'a'][1] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'a'][2] = gPHMM_syn;
	gPHMM_ALIGN_SCORES[(int)'a'][3] = gPHMM_nsyn;

	gPHMM_ALIGN_SCORES[(int)'c'][0] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'c'][1] = gPHMM_match;
	gPHMM_ALIGN_SCORES[(int)'c'][2] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'c'][3] = gPHMM_syn;

	gPHMM_ALIGN_SCORES[(int)'g'][0] = gPHMM_syn;
	gPHMM_ALIGN_SCORES[(int)'g'][1] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'g'][2] = gPHMM_match;
	gPHMM_ALIGN_SCORES[(int)'g'][3] = gPHMM_nsyn;

	gPHMM_ALIGN_SCORES[(int)'t'][0] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'t'][1] = gPHMM_syn;
	gPHMM_ALIGN_SCORES[(int)'t'][2] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'t'][3] = gPHMM_match;

	gPHMM_ALIGN_SCORES[(int)'n'][0] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'n'][1] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'n'][2] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'n'][3] = gPHMM_nsyn;
	//*/

}


/*
 * Set them back to the way they were at the beginning
 */
void reset_alignment_matrices() {
	// If we've specified something, return without resetting
	if(pos_matrix)
		return;
		
	gALIGN_SCORES[(int)'a'][0] = gALIGN_SCORES[(int)'A'][0] = gMATCH;
	gALIGN_SCORES[(int)'a'][1] = gALIGN_SCORES[(int)'A'][1] = gTRANSVERSION;
	gALIGN_SCORES[(int)'a'][2] = gALIGN_SCORES[(int)'A'][2] = gTRANSITION;
	gALIGN_SCORES[(int)'a'][3] = gALIGN_SCORES[(int)'A'][3] = gTRANSVERSION;

	gALIGN_SCORES[(int)'c'][0] = gALIGN_SCORES[(int)'C'][0] = gTRANSVERSION;
	gALIGN_SCORES[(int)'c'][1] = gALIGN_SCORES[(int)'C'][1] = gMATCH;
	gALIGN_SCORES[(int)'c'][2] = gALIGN_SCORES[(int)'C'][2] = gTRANSVERSION;
	gALIGN_SCORES[(int)'c'][3] = gALIGN_SCORES[(int)'C'][3] = gTRANSITION;

	gALIGN_SCORES[(int)'g'][0] = gALIGN_SCORES[(int)'G'][0] = gTRANSITION;
	gALIGN_SCORES[(int)'g'][1] = gALIGN_SCORES[(int)'G'][1] = gTRANSVERSION;
	gALIGN_SCORES[(int)'g'][2] = gALIGN_SCORES[(int)'G'][2] = gMATCH;
	gALIGN_SCORES[(int)'g'][3] = gALIGN_SCORES[(int)'G'][3] = gTRANSVERSION;

	gALIGN_SCORES[(int)'t'][0] = gALIGN_SCORES[(int)'T'][0] = gTRANSVERSION;
	gALIGN_SCORES[(int)'t'][1] = gALIGN_SCORES[(int)'T'][1] = gTRANSITION;
	gALIGN_SCORES[(int)'t'][2] = gALIGN_SCORES[(int)'T'][2] = gTRANSVERSION;
	gALIGN_SCORES[(int)'t'][3] = gALIGN_SCORES[(int)'T'][3] = gMATCH;
	

	
	gPHMM_ALIGN_SCORES[(int)'a'][0] = gPHMM_ALIGN_SCORES[(int)'A'][0] = gPHMM_match;
	gPHMM_ALIGN_SCORES[(int)'a'][1] = gPHMM_ALIGN_SCORES[(int)'A'][1] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'a'][2] = gPHMM_ALIGN_SCORES[(int)'A'][2] = gPHMM_syn;
	gPHMM_ALIGN_SCORES[(int)'a'][3] = gPHMM_ALIGN_SCORES[(int)'A'][3] = gPHMM_nsyn;

	gPHMM_ALIGN_SCORES[(int)'c'][0] = gPHMM_ALIGN_SCORES[(int)'C'][0] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'c'][1] = gPHMM_ALIGN_SCORES[(int)'C'][1] = gPHMM_match;
	gPHMM_ALIGN_SCORES[(int)'c'][2] = gPHMM_ALIGN_SCORES[(int)'C'][2] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'c'][3] = gPHMM_ALIGN_SCORES[(int)'C'][3] = gPHMM_syn;

	gPHMM_ALIGN_SCORES[(int)'g'][0] = gPHMM_ALIGN_SCORES[(int)'G'][0] = gPHMM_syn;
	gPHMM_ALIGN_SCORES[(int)'g'][1] = gPHMM_ALIGN_SCORES[(int)'G'][1] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'g'][2] = gPHMM_ALIGN_SCORES[(int)'G'][2] = gPHMM_match;
	gPHMM_ALIGN_SCORES[(int)'g'][3] = gPHMM_ALIGN_SCORES[(int)'G'][3] = gPHMM_nsyn;

	gPHMM_ALIGN_SCORES[(int)'t'][0] = gPHMM_ALIGN_SCORES[(int)'T'][0] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'t'][1] = gPHMM_ALIGN_SCORES[(int)'T'][1] = gPHMM_syn;
	gPHMM_ALIGN_SCORES[(int)'t'][2] = gPHMM_ALIGN_SCORES[(int)'T'][2] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'t'][3] = gPHMM_ALIGN_SCORES[(int)'T'][3] = gPHMM_match;

	gPHMM_ALIGN_SCORES[(int)'n'][0] = gPHMM_ALIGN_SCORES[(int)'N'][0] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'n'][1] = gPHMM_ALIGN_SCORES[(int)'N'][1] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'n'][2] = gPHMM_ALIGN_SCORES[(int)'N'][2] = gPHMM_nsyn;
	gPHMM_ALIGN_SCORES[(int)'n'][3] = gPHMM_ALIGN_SCORES[(int)'N'][3] = gPHMM_nsyn;
}
