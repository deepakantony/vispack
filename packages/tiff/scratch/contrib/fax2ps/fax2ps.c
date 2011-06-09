#ifndef lint
static char rcsid[] = "$Header: /uusoc/res/image/CVS/vispack/packages/tiff/scratch/contrib/fax2ps/fax2ps.c,v 1.1.1.1 2003-02-12 16:51:50 whitaker Exp $";
#endif

/*
 * Copyright (c) 1991, 1992 by Sam Leffler.
 * All rights reserved.
 *
 * This file is provided for unrestricted use provided that this
 * legend is included on all tape media and as a part of the
 * software program in whole or part.  Users may copy, modify or
 * distribute this file at will.
 */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "defs.h"

#define	MAXCODEPROBES	5
CodeEntry* codehash[CODEHASH];
CodeEntry codetable[MAXCODES];

char**	codeNames;		/* codeNames[code] => ASCII code name */
char*	codeNameSpace;		/* storage space for codeNames and strings */
int	ncodes = 0;		/* number of assigned codes */
int	includeStatistics = FALSE;/* if 1, add comments w/ frequency stats */
int	startOfRow;		/* if 1, have yet to emit a code for this row */
int	dopairs = FALSE;	/* if 1, encode pairs of codes */
float	defxres = 204.;		/* default x resolution (pixels/inch) */
float	defyres = 98.;		/* default y resolution (lines/inch) */

CodeEntry*
DECLARE2(enterCode, int, dx, int, len)
{
    int h, c;
    CodeEntry* cp;

    cp = codehash[h = HASHCODE(dx,len)];
    if (cp) {
	if (cp->move == dx && cp->runlen == len)
	    return (cp);
	c = h ? CODEHASH - h : 1;	/* Knott's rehash algorithm */
	for (;;) {
	    if ((h -= c) < 0)
		h += CODEHASH;
	    cp = codehash[h];
	    if (!cp)
		break;
	    if (cp->move == dx && cp->runlen == len)
		return (cp);
	}
    }
    if (ncodes == MAXCODES) {
	fprintf(stderr, "Panic, code table overflow\n");
	exit(-1);
    }
    codehash[h] = cp = &codetable[ncodes++];
    cp->c.count = 0;
    cp->c.code = (u_short) -1;
    cp->c.cost = 0;
    cp->move = dx;
    cp->runlen = len;
    return (cp);
}

void
DECLARE2(printCode, TIFF*, tif, CodeEntry*, cp)
{
    if (startOfRow) {
	printf("%d r\n", fax.row);
	startOfRow = FALSE;
    }
    if (cp->c.code == (u_short) -1)
	printf("%d %d f", cp->runlen, cp->move);
    else
	printf("%s", codeNames[cp->c.code]);
    putchar('\n');
}

CodePairEntry* pairhash[PAIRHASH];
CodePairEntry pairtable[MAXPAIRS];
int	npairs = 0;

CodePairEntry*
DECLARE2(findPair, CodeEntry*, a, CodeEntry*, b)
{
    int h, c;
    CodePairEntry* pp;

    pp = pairhash[h = HASHPAIR(a,b)];
    if (pp->a == a && pp->b == b)
	return (pp);
    c = h ? PAIRHASH - h : 1;		/* Knott's rehash algorithm */
    for (;;) {
	if ((h -= c) < 0)
	    h += PAIRHASH;
	pp = pairhash[h];
	if (pp->a == a && pp->b == b)
	    return (pp);
    }
    /*NOTREACHED*/
}

CodePairEntry*
DECLARE2(enterPair, CodeEntry*, a, CodeEntry*, b)
{
    int h, c;
    CodePairEntry* pp;

    pp = pairhash[h = HASHPAIR(a,b)];
    if (pp) {
	if (pp->a == a && pp->b == b)
	    return (pp);
	c = h ? PAIRHASH - h : 1;	/* Knott's rehash algorithm */
	for (;;) {
	    if ((h -= c) < 0)
		h += PAIRHASH;
	    pp = pairhash[h];
	    if (!pp)
		break;
	    if (pp->a == a && pp->b == b)
		return (pp);
	}
    }
    if (npairs == MAXPAIRS) {
	fprintf(stderr, "Help, pair table overflow\n");
	exit(-1);
    }
    pairhash[h] = pp = &pairtable[npairs++];
    pp->c.count = 0;
    pp->c.code = (u_short) -1;
    pp->c.cost = 0;
    pp->a = a;
    pp->b = b;
    return (pp);
}

int
DECLARE3(printPair, TIFF*, tif, CodeEntry*, a, CodeEntry*, b)
{
    CodePairEntry* pp = findPair(a,b);

    if (pp && pp->c.code != (u_short) -1) {
	if (startOfRow) {
	    printf("%d r\n", fax.row);
	    startOfRow = FALSE;
	}
	printf("%s\n", codeNames[pp->c.code]);
	return (TRUE);
    } else
	return (FALSE);
}

#define	MIN(a,b)	((a)<(b)?(a):(b))

static char alphabet[] = {
    '!', '"', '#', '$', '&', '*', '+', ',', ':', ';', '?', '@',
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
    'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
    '^', '_', '~'
};
static int MaxAlpha = (sizeof (alphabet) / sizeof (alphabet[0]));

/*
 * Construct a table of code names from a limited
 * alphabet.  By default, the alphabet is comprised
 * of all the non-special ASCII PostScript characters
 * (except for lower case alphabetics which are avoided
 * to insure we don't construct a PostScript operator
 * name).  The alphabet can be restricted to just the 
 * upper case alphabetics with the -a option.
 */
void
DECLARE1(makeCodeNames, int, ncodes)
{
    int n, cc, len, nc, digit, code;
    char* cp;
    short key[11];

    n = ncodes;
    cc = 0;
    len = 2;		/* minimum length string */
    nc = 1;
    do {
	nc *= MaxAlpha;
	cc += MIN(n, nc) * len++;
    } while ((n -= MIN(n,nc)) > 0);
    codeNameSpace = (char*) malloc(ncodes*sizeof (char*) + cc);
    codeNames = (char**) codeNameSpace;
    cp = (char*) (codeNames + ncodes);

    n = ncodes;
    len = 1;
    code = 0;
    nc = 1;
    for (digit = 0; digit < 11; digit++)
	key[digit] = -1;
    do {
	nc *= MaxAlpha;
	cc = MIN(n,nc);
	while (cc-- > 0) {
	    for (digit = 0; ++key[digit] == MaxAlpha; key[digit++] = 0)
		;
	    codeNames[code++] = cp;
	    switch (len) {
	    case 5: *cp++ = alphabet[key[4]];
	    case 4: *cp++ = alphabet[key[3]];
	    case 3: *cp++ = alphabet[key[2]];
	    case 2: *cp++ = alphabet[key[1]];
	    case 1: *cp++ = alphabet[key[0]];
	    }
	    *cp++ = '\0';
	}
	len++;
    } while ((n -= MIN(n,nc)) > 0);
}

#define	DIGITS(x)	((x) < 10 ? 1 : (x) < 100 ? 2 : (x) < 1000 ? 3 : 4)

int
DECLARE1(codeCost, CodeEntry*, cp)
{
    /* 3 is constant overhead for <blank><len><blank><move>f */
    return (3 + DIGITS(cp->move) + DIGITS(cp->runlen));
}

/*
 * DEFCOST is the cost to define the specified code for this
 * move-draw operation.  It's calculated according to:
 *	/<code>{<move-draw>}d\n
 * which translates to 5 + cost(code) + cost(move-draw).
 */
#define	DEFCOST(mdc, codelen)	(mdc + codelen + 5)
/*
 * USECOST is the cost to use the code once it's been defined.
 * This is just the number of uses times the length of the code
 * plus one -- the one is for the blank needed by the parser.
 */
#define	USECOST(count, codelen)	(count * (codelen+1))
/*
 * CODEDCOST is the total cost to define and use the defined
 * code (in place of straight move-draw operations).  We assume
 * that the "code" field holds the result of the {pair|code}Cost
 * calculation (see below).
 */
#define	CODEDCOST(c, codelen) \
    (USECOST(c->count, codelen) + DEFCOST(c->code, codelen))

int
DECLARE1(codeLength, int, code)
{
    return (code < MaxAlpha ? 1 : code < MaxAlpha*MaxAlpha ? 2 : 3);
}

int
DECLARE2(codeCostCompare, Code**, c1, Code**, c2)
{
    return ((*c2)->cost - (*c1)->cost);
}

/*
 * Sort pairs and singletons, assign codes, and 
 * reset counts for the next pass where we test
 * out the encoding.
 */
int
DECLARE2(assignCodes, Code**, sorted, int, clearCounts)
{
    int i, n, code, codelen;
    CodePairEntry* pp;
    CodeEntry* cp;

    /*
     * Calculate the cost to use the operation "as is"
     * (i.e. w/o an encoding).  This is simply the number
     * of uses times the cost of the move-draw operation.
     * We save the {pair|code}Cost calculation for reuse
     * below (when making the final decision about which
     * codes are actually worth defining.
     */
    i = 0;
    for (cp = codetable, n = ncodes; n-- > 0; i++, cp++) {
	sorted[i] = &cp->c;
	cp->c.code = codeCost(cp);
	cp->c.cost = cp->c.count * cp->c.code;
    }
    for (pp = pairtable, n = npairs; n-- > 0; i++, pp++) {
	sorted[i] = &pp->c;
	pp->c.code = pp->a->c.code + pp->b->c.code;
	pp->c.cost = pp->c.count * pp->c.code;
    }
    qsort(sorted, ncodes + npairs, sizeof (Code*), codeCostCompare);
    code = 0;
    codelen = codeLength(code);
    for (i = 0, n = npairs + ncodes; n-- > 0; i++) {
	Code* c = sorted[i];
	if (c->cost > CODEDCOST(c, codelen)) {
	    c->code = code++;
	    codelen = codeLength(code);
	} else
	    c->code = (u_short) -1;
	if (clearCounts)
	    c->count = 0;
    }
    return (code);
}

void
makeCodeTable()
{
    int i, n;
    Code** sorted;
    CodeEntry* cp;

    sorted = (Code**) malloc((npairs+ncodes) * sizeof (Code*));
    makeCodeNames(n = assignCodes(sorted, FALSE));
    /*
     * Oversize the dictionary in case the lookup function
     * is optimized for a partially populated data structure.
     * Most any hashing algorithm should do fine with a max
     * 75% population.
     */
    printf("%d dict begin\n", 4*n/3);
    for (i = 0, n = npairs + ncodes; n-- > 0; i++) {
	Code* c = sorted[i];
	if (c->code == (u_short) -1)
	    continue;
	if (isPair(c)) {
	    CodePairEntry* pp = (CodePairEntry*) c;
	    printf("/%s{", codeNames[c->code]);
	    cp = pp->a;
	    if (cp->c.code < c->code)
		printf("%s", codeNames[cp->c.code]);
	    else
		printf("%d %d f", cp->runlen, cp->move);
	    putchar(' ');
	    cp = pp->b;
	    if (cp->c.code < c->code)
		printf("%s", codeNames[cp->c.code]);
	    else
		printf("%d %d f", cp->runlen, cp->move);
	    printf("}d");
	    if (includeStatistics)
		printf("\t%% %d hits", c->count);
	    putchar('\n');
	} else {
	    cp = (CodeEntry*) c;
	    printf("/%s{%d %d f}d", codeNames[c->code], cp->runlen, cp->move);
	    if (includeStatistics)
		printf("\t%% %d hits", c->count);
	    putchar('\n');
	}
    }
    free((char*) sorted);
}

static void
DECLARE2(setupPass, TIFF*, tif, int, pass)
{
    long* stripbytecount;

    fax.pass = pass;
    fax.row = 0;
    TIFFGetField(tif, TIFFTAG_STRIPBYTECOUNTS, &stripbytecount);
    fax.cc = stripbytecount[0];
    fax.bp = fax.buf;
    FaxPreDecode(tif);
}

static	int totalPages = 0;

void
DECLARE2(printTIF, TIFF*, tif, int, pageNumber)
{
    u_long w, h;
    short fill, unit, photometric, compression;
    float xres, yres;
    long g3opts;
    long* stripbytecount;

    TIFFGetField(tif, TIFFTAG_STRIPBYTECOUNTS, &stripbytecount);
    fax.cc = stripbytecount[0];
    fax.buf = (u_char*) malloc(fax.cc);
    TIFFReadRawStrip(tif, 0, fax.buf, fax.cc);
    if (TIFFGetField(tif,TIFFTAG_FILLORDER, &fill) && fill != FILLORDER_MSB2LSB)
	TIFFReverseBits(fax.buf, fax.cc);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
    if (!TIFFGetField(tif, TIFFTAG_XRESOLUTION, &xres)) {
	TIFFWarning(TIFFFileName(tif),
	    "No x-resolution, assuming %g dpi", defxres);
	xres = defxres;
    }
    if (!TIFFGetField(tif, TIFFTAG_YRESOLUTION, &yres)) {
	TIFFWarning(TIFFFileName(tif),
	    "No y-resolution, assuming %g lpi", defyres);
	yres = defyres;					/* XXX */
    }
    if (TIFFGetField(tif, TIFFTAG_RESOLUTIONUNIT, &unit) &&
      unit == RESUNIT_CENTIMETER) {
	xres *= 25.4;
	yres *= 25.4;
    }
    TIFFGetField(tif, TIFFTAG_PHOTOMETRIC, &photometric);
    fax.b.white = (photometric == PHOTOMETRIC_MINISBLACK);
    /*
     * Calculate the scanline/tile widths.
     */
    if (isTiled(tif)) {
        fax.b.rowbytes = TIFFTileRowSize(tif);
	TIFFGetField(tif, TIFFTAG_TILEWIDTH, &fax.b.rowpixels);
    } else {
        fax.b.rowbytes = TIFFScanlineSize(tif);
	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &fax.b.rowpixels);
    }
    if (TIFFGetField(tif, TIFFTAG_GROUP3OPTIONS, &g3opts))
	fax.is2d = (g3opts & GROUP3OPT_2DENCODING) != 0;
    else
	fax.is2d = 0;
    fax.scanline = (u_char*) malloc(2*fax.b.rowbytes);
    fax.b.refline = fax.scanline + fax.b.rowbytes;
    TIFFGetField(tif, TIFFTAG_COMPRESSION, &compression);
    if (compression == COMPRESSION_CCITTFAX4)
	fax.options = FAX3_NOEOL;
    /*
     * First pass: create code and pair frequency tables.
     */
    bzero((char*) codehash, sizeof (codehash));
    bzero((char*) pairhash, sizeof (pairhash));
    setupPass(tif, 1);
    ncodes = npairs = 0;
    if (compression == COMPRESSION_CCITTFAX4)
	Fax4DecodeRow(tif, h*w);
    else
	Fax3DecodeRow(tif, h*w);
    if (dopairs) {
	/*
	 * Second pass: assign codes according to frequency
	 *  thresholds and rescan data to get true usage.
	 */
	CodeEntry** sorted = (CodeEntry**)
	    malloc((npairs+ncodes) * sizeof (CodeEntry*));
	assignCodes((Code**) sorted, TRUE);
	free((char*) sorted);
	setupPass(tif, 2);
	if (compression == COMPRESSION_CCITTFAX4)
	    Fax4DecodeRow(tif, h*w);
	else
	    Fax3DecodeRow(tif, h*w);
    }
    /*
     * Third pass: generate code table and encoded data.
     */
    printf("%%%%Page: \"%d\" %d\n", pageNumber, pageNumber);
    printf("gsave\n");
    printf("0 %d translate\n", (int)(h/yres*72.));
    printf("%g %g scale\n", 72./xres, -72./yres);
    printf("0 setgray\n");
    makeCodeTable();
    setupPass(tif, 3);
    while (fax.cc > 0) {
	startOfRow = 1;
	if (compression == COMPRESSION_CCITTFAX4)
	    Fax4DecodeRow(tif, w);
	else
	    Fax3DecodeRow(tif, w);
	if (!startOfRow)
	    printf("s\n");
    }
    printf("p\n");
    printf("grestore\n");
    free((char*) fax.scanline);
    free((char*) fax.buf);
    free(codeNameSpace);
    totalPages++;
}

#define	GetPageNumber(tif) \
TIFFGetField(tif, TIFFTAG_PAGENUMBER, &pn, &ptotal)

int
DECLARE2(findPage, TIFF*, tif, int, pageNumber)
{
    short pn = -1, ptotal = -1;
    if (GetPageNumber(tif)) {
	while (pn != pageNumber && TIFFReadDirectory(tif) && GetPageNumber(tif))
	    ;
	return (pn == pageNumber);
    } else
	return (TIFFSetDirectory(tif, pageNumber-1));
}

void
DECLARE4(fax2ps, TIFF*, tif, int, npages, int*, pages, char*, filename)
{
    if (npages > 0) {
	short pn, ptotal;
	int i;

	if (!GetPageNumber(tif))
	    fprintf(stderr, "%s: No page numbers, counting directories.\n",
		filename);
	for (i = 0; i < npages; i++) {
	    if (findPage(tif, pages[i]))
		printTIF(tif, pages[i]);
	    else
		fprintf(stderr, "%s: No page number %d\n", filename, pages[i]);
	}
    } else {
	int pageNumber = 1;
	do
	    printTIF(tif, pageNumber++);
	while (TIFFReadDirectory(tif));
    }
}

#undef GetPageNumber

static int
DECLARE2(pcompar, void*, va, void*, vb)
{
    int* pa = (int*) va;
    int* pb = (int*) vb;
    return (*pa - *pb);
}

extern	double atof();

DECLARE2(main, int, argc, char**, argv)
{
    extern int optind;
    extern char* optarg;
    int c, pageNumber;
    int* pages = 0, npages = 0;
    int dowarnings = FALSE;	/* if 1, enable library warnings */
    long t;
    TIFF* tif;

    while ((c = getopt(argc, argv, "p:x:y:aswz")) != -1)
	switch (c) {
	case 'a':		/* use only upper-case alphabetics */
	    MaxAlpha = 26;
	    for (t = 0; t < MaxAlpha; t++)
		alphabet[t] = 'A' + t;
	    break;
	case 'p':		/* print specific page */
	    pageNumber = atoi(optarg);
	    if (pageNumber < 1) {
		fprintf(stderr, "%s: Invalid page number (must be > 0).\n",
		    optarg);
		exit(-1);
	    }
	    if (pages)
		pages = (int*) realloc((char*) pages, (npages+1)*sizeof (int));
	    else
		pages = (int*) malloc(sizeof (int));
	    pages[npages++] = pageNumber;
	    break;
	case 's':		/* include frequency statistics as comments */
	    includeStatistics = TRUE;
	    break;
	case 'w':
	    dowarnings = TRUE;
	    break;
	case 'x':
	    defxres = atof(optarg);
	    break;
	case 'y':
	    defyres = atof(optarg);
	    break;
	case 'z':		/* do pair analysis (not effective) */
	    dopairs = TRUE;
	    break;
	case '?':
	    fprintf(stderr,
"usage: %s [-a] [-w] [-p pagenumber] [-x xres] [-y res] [-s] [files]\n",
		argv[0]);
	    exit(-1);
	}
    if (npages > 0)
	qsort(pages, npages, sizeof (int), pcompar);
    if (!dowarnings)
	TIFFSetWarningHandler(0);
    printf("%%!PS-Adobe-3.0\n");
    printf("%%%%Creator: fax2ps\n");
#ifdef notdef
    printf("%%%%Title: %s\n", file);
#endif
    t = time(0);
    printf("%%%%CreationDate: %s", ctime(&t));
    printf("%%%%Origin: 0 0\n");
    printf("%%%%BoundingBox: 0 0 %d %d\n", 11*72, (int)(8.5*72));/* XXX */
    printf("%%%%Pages: (atend)\n");
    printf("%%%%EndComments\n");
    printf("%%%%BeginProlog\n");
    printf("/d{bind def}def\n");			/* bind and def proc */
    printf("/f{0 rmoveto 0 rlineto}d\n");		/* fill span */
    printf("/r{0 exch moveto}d\n");			/* begin row */
    printf("/s{stroke}d\n");				/* stroke row */
    printf("/p{end showpage}d \n");			/* end page */
    printf("%%%%EndProlog\n");
    if (optind < argc) {
	for (; optind < argc; optind++) {
	    tif = TIFFOpen(argv[optind], "r");
	    if (!tif) {
		fprintf(stderr, "%s: Can not open, or not a TIFF file.\n",
		    argv[optind]);
		continue;
	    }
	    fax2ps(tif, npages, pages, argv[optind]);
	    TIFFClose(tif);
	}
    } else {
	int n, fd;
	char temp[1024], buf[16*1024];

	strcmp(temp, "/tmp/fax2psXXXXXX");
	fd = mkstemp(temp);
	if (fd == -1) {
	    fprintf(stderr, "Could not create temp file \"%s\"\n", temp);
	    exit(-2);
	}
	while ((n = read(fileno(stdin), buf, sizeof (buf))) > 0)
	    write(fd, buf, n);
	tif = TIFFOpen(temp, "r");
	unlink(temp);
	if (tif) {
	    fax2ps(tif, npages, pages, "<stdin>");
	    TIFFClose(tif);
	} else
	    fprintf(stderr, "%s: Can not open, or not a TIFF file.\n", temp);
	close(fd);
    }
    printf("%%%%Trailer\n");
    printf("%%%%Pages: %u\n", totalPages);
    printf("%%%%EOF\n");
    exit(0);
}
