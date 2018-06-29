
/*
 *  vtprogressbar.c
 *  implementation for a very simple
 *  progress bar
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/05/06 02:11:30 CET
 *  
 *  SVN
 *  Revision of last commit: $Rev: 19 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-05-14 15:43:29 +0200 (Wed, 14 May 2008) $
 *
 *  Id: $Id: vtprogressbar.c 19 2008-05-14 13:43:29Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/branches/esa/trunk/libs/vtprogressbar.c $
 */
 #include <stdio.h>
 #include <stdlib.h>
 #include "basic-types.h"
 #include "vtprogressbar.h"


void
cursorInvisible() {
  fprintf(stderr, "%c%c%c%d%c", 27, '[', '?', 25, 'l');
}

void
cursorVisible() {
  fprintf(stderr, "%c%c%c%d%c", 27, '[', '?', 25, 'h');
}



/*---------------------------- initProgressBarVT -----------------------------
 *    
 * initializes a progress bar for VT
 * 
 */
 
void
initProgressBarVT ()
{	
  	fprintf(stderr, "%c%c%c", 27, '[', 's');
    fprintf(stderr, "%c%c%c", 27, '[', 'K');
    return ;
}

/*------------------------------ progressBarVT -------------------------------
 *    
 * a simple progress bar for VT terminals
 * 
 */

void
progressBarVT (char *message, Uint complete, Uint processed, Uint size)
{
    Uint i, percent, bar;
	char cur;
    
    if (complete == 0) complete = 1;
	bar=(size*processed)/(complete);
	percent =(processed*100)/(complete);
	fprintf(stderr, "[");
	for(i=0; i < size; i++) {
		if(i<=bar) fprintf(stderr, "=");
		else fprintf(stderr," ");
	}
	i = processed % 30;
	if (i<=10) cur = '/'; 
	else if (i<=20)cur='\\';
	else cur='-';
	fprintf(stderr,"]   %d%c(%d)  %s  %c\n", percent, '%', processed, message, cur);
	fprintf(stderr,"%c%c%c", 27, '[', 'A');
	return;
}


