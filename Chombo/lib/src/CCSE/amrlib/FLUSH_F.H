#ifndef _FLUSH_F_H_
#define _FLUSH_F_H_
c
c $Id: FLUSH_F.H,v 1.1 2007-05-24 23:52:09 tdsternberg Exp $
c

#if defined(BL_AIX)
#define FORT_FLUSH flush_
#else
#define FORT_FLUSH flush
#endif

#endif

