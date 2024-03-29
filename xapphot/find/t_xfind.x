include <fset.h>
include <imhdr.h>
include "../lib/xphot.h"
include "../lib/display.h"
include "../lib/objects.h"


# T_XFIND -- Interactively display a list of images.

procedure t_xfind()

int	dirlist, imlist, objlist, reslist, greslist, imno, olno, rlno, glno
int	ol, rl, gl, wcs, key, nobjs, update, verbose, iupdate
pointer	sp, images, objects, results, graphics, imname, olname, rlname, glname
pointer	cmd, str, pstatus
pointer	gd, xp, im, symbol, lsymbol
real	owx, owy, wx, wy
bool	clgetb()
int	xp_dirlist(), xp_imlist(), imtrgetim(), imtlen(), xp_stati(), clgcur()
int	btoi(), fntlenb(), xp_mkolist(), fntrfnb(), open(), fstati()
int	xp_mkrlist(), access(), xp_afind(), xp_pmodfit()
pointer	gopen(), immap(), xp_dmeas(), stfind, stenter(), xp_statp()
pointer	xp_xdcolon()
errchk	immap(), open()

define	noninteractive_ 99

begin
	call fseti (STDOUT, F_FLUSHNL, YES)

	call smark (sp)
	call salloc (images, SZ_FNAME, TY_CHAR)
	call salloc (objects, SZ_FNAME, TY_CHAR)
	call salloc (results, SZ_FNAME, TY_CHAR)
	call salloc (graphics, SZ_FNAME, TY_CHAR)
	call salloc (imname, SZ_FNAME, TY_CHAR)
	call salloc (olname, SZ_FNAME, TY_CHAR)
	call salloc (rlname, SZ_FNAME, TY_CHAR)
	call salloc (glname, SZ_FNAME, TY_CHAR)
	call salloc (cmd, SZ_LINE, TY_CHAR)
	call salloc (str, SZ_FNAME, TY_CHAR)

	# Get the task parameters.
	call clgstr ("images", Memc[images], SZ_FNAME)
	call clgstr ("objects", Memc[objects], SZ_FNAME)
	call clgstr ("results", Memc[results], SZ_FNAME)
	call clgstr ("graphics", Memc[graphics], SZ_FNAME)

	# Initialize the algorithm parameters.
	call xp_gxdpars (xp)
	pstatus = xp_statp(xp, PSTATUS)

	# Open the graphics stream.
	if (clgetb ("interactive")) {
	    if (xp_stati (xp, DERASE) == YES)
	        gd = gopen (Memc[graphics], NEW_FILE, STDGRAPH)
	    else
	        gd = gopen (Memc[graphics], APPEND, STDGRAPH)
	} else
	    gd = NULL
	update = btoi (clgetb ("update"))
	LOGRESULTS(pstatus) = btoi (clgetb("logresults"))
	verbose = btoi (clgetb ("verbose"))

	# Open the current directory list.
	dirlist = xp_dirlist ("..,*")
	call fpathname ("", Memc[cmd], SZ_LINE)
	call xp_sets (xp, STARTDIR, Memc[cmd])
	call xp_sets (xp, CURDIR, Memc[cmd])

	# Open the image list.
	imlist = xp_imlist (Memc[images])
	call xp_sets (xp, IMTEMPLATE, Memc[images])
	imno = 1
	NEWIMAGE(pstatus) = YES

	# Open the objects lists.
	objlist = xp_mkolist (imlist, Memc[objects], "default", "obj")
	call xp_sets (xp, OFTEMPLATE, Memc[objects])
	olno = 1
	NEWLIST(pstatus) = YES

	# Open the results file list.
	reslist = xp_mkrlist (imlist, Memc[results], "default", "fnd", NO)
	call xp_sets (xp, RFTEMPLATE, Memc[results])
	rlno = 1
	call xp_oseqlist (xp)
	greslist = xp_mkrlist (imlist, "", "default", "geo", NO)
	call xp_sets (xp, GFTEMPLATE, "")
	glno = 1
	NEWRESULTS(pstatus) = YES

	repeat {

	    # Open the image and display it.
	    if (NEWIMAGE(pstatus) == YES) {
	        if (imtrgetim (imlist, imno, Memc[imname], SZ_FNAME) == EOF) {
		    im = NULL
		    Memc[imname] = EOS
	            call xp_seti (xp, IMNUMBER, 0)
	        } else { 
	            iferr {
	                im = immap (Memc[imname], READ_ONLY, 0)
	            } then {
			if (gd != NULL) {
		            call gclear (gd)
		            call gflush (gd)
			}
		        #call printf (
			#"Warning: Cannot open image (%s)\n")
			    #call pargstr (Memc[imname])
		        im = NULL
	            } else {
			if (gd != NULL)
	                    call xp_display (gd, xp, im, 1, IM_LEN(im,1), 1,
		                IM_LEN(im,2), IMAGE_DISPLAY_WCS,
				IMAGE_DISPLAY_WCS)
	            }
		    call xp_keyset (im, xp)
	            call xp_seti (xp, IMNUMBER, imno)
	        }
	        call xp_sets (xp, IMAGE, Memc[imname])
	    }

	    # Open the coordinate file if any.
	    if (NEWLIST(pstatus) == YES) {
		OBJNO(pstatus) = 0
		call xp_clsobjects (xp)
	        if (im == NULL || fntlenb (objlist) <= 0) {
		    Memc[olname] = EOS
	            call xp_seti (xp, OFNUMBER, 0)
		} else if (NEWIMAGE(pstatus) == NO && fntrfnb (objlist, olno,
		    Memc[olname], SZ_FNAME) != EOF) {
	            call xp_seti (xp, OFNUMBER, olno)
	        } else if (fntrfnb (objlist, xp_stati(xp,IMNUMBER),
		    Memc[olname], SZ_FNAME) != EOF) {
	            call xp_seti (xp, OFNUMBER, xp_stati (xp, IMNUMBER))
	        } else if (fntrfnb (objlist, fntlenb (objlist),
		    Memc[olname], SZ_FNAME) != EOF) {
	            call xp_seti (xp, OFNUMBER, fntlenb (objlist))
		}
		if (Memc[olname] == EOS) {
		    ol = NULL
	        } else iferr {
	            ol = open (Memc[olname], READ_ONLY, TEXT_FILE)
	        } then {
		    #call printf (
		    #"Warning: Cannot open objects file (%s)\n")
		        #call pargstr (Memc[olname])
		    ol = NULL
		}
	        call xp_sets (xp, OBJECTS, Memc[olname])
	    }

	    # Open the results file if any.
	    if (NEWRESULTS(pstatus) == YES) {
	        if (im == NULL || fntlenb (reslist) <= 0) {
		    Memc[rlname] = EOS
	            call xp_seti (xp, RFNUMBER, 0)
		} else if (NEWIMAGE(pstatus) == NO && fntrfnb (reslist, rlno,
		    Memc[rlname], SZ_FNAME) != EOF) {
	            call xp_seti (xp, RFNUMBER, rlno)
	        } else if (fntrfnb (reslist, xp_stati(xp,IMNUMBER),
		    Memc[rlname], SZ_FNAME) != EOF) {
	            call xp_seti (xp, RFNUMBER, xp_stati (xp, IMNUMBER))
	        } else if (fntrfnb (reslist, fntlenb (reslist),
		    Memc[rlname], SZ_FNAME) != EOF) {
	            call xp_seti (xp, RFNUMBER, fntlenb (reslist))
		}
		if (Memc[rlname] == EOS) {
		    rl = NULL
	        } else iferr {
		    if (access (Memc[rlname], 0, 0) == YES) {
	                rl = open (Memc[rlname], APPEND, TEXT_FILE)
                        lsymbol = stfind (xp_statp(xp, SEQNOLIST), Memc[rlname])
                        if (lsymbol == NULL)
                            SEQNO(pstatus) = 0
                        else
                            SEQNO(pstatus) = XP_MAXSEQNO(lsymbol)
		    } else {
	                rl = open (Memc[rlname], NEW_FILE, TEXT_FILE)
			SEQNO(pstatus) = 0
		    }
		    if (SEQNO(pstatus) == 0)
			call xp_whfind (xp, rl, "xfind")
		    call xp_whiminfo (xp, rl)
		    if (SEQNO(pstatus) == 0)
			call xp_xfbnr (xp, rl)
	        } then {
		    #call printf (
		    #"Warning: Cannot open objects file (%s)\n")
		        #call pargstr (Memc[rlname])
		    rl = NULL
		}
	        call xp_sets (xp, RESULTS, Memc[rlname])

                if (im == NULL || fntlenb (greslist) <= 0) {
                    Memc[glname] = EOS
                    call xp_seti (xp, GFNUMBER, 0)
                } else if (NEWIMAGE(pstatus) == NO && fntrfnb (greslist, glno,
                    Memc[glname], SZ_FNAME) != EOF) {
                    call xp_seti (xp, GFNUMBER, glno)
                } else if (fntrfnb (greslist, xp_stati(xp,IMNUMBER),
                    Memc[glname], SZ_FNAME) != EOF) {
                    call xp_seti (xp, GFNUMBER, xp_stati (xp, IMNUMBER))
                } else if (fntrfnb (greslist, fntlenb (greslist),
                    Memc[glname], SZ_FNAME) != EOF) {
                    call xp_seti (xp, GFNUMBER, fntlenb (greslist))
                }
                if (Memc[glname] == EOS) {
                    gl = NULL
                } else iferr {
                    if (access (Memc[glname], 0, 0) == YES)
                        gl = open (Memc[glname], APPEND, TEXT_FILE)
                    else
                        gl = open (Memc[glname], NEW_FILE, TEXT_FILE)
                } then {
                    #call printf (
                    #"Warning: Cannot open objects file (%s)\n")
                        #call pargstr (Memc[glname])
                    gl = NULL
                }
                call xp_sets (xp, GRESULTS, Memc[glname])
	    }


	    if (gd == NULL)
		goto noninteractive_

	    # Initialize the analysis.
	    NEWIMAGE(pstatus) = NO; NEWLIST(pstatus) = NO;
	    NEWRESULTS(pstatus) = NO
	    symbol = NULL; owx = INDEFR; owy = INDEFR

	    # Loop over the cursor commands.
	    while (clgcur ("gcommands", wx, wy, wcs, key, Memc[cmd],
		SZ_LINE) != EOF) {


		switch (key) {

		# Quit the program.
		case 'Q':
		    imno = EOF
		    break

		# Erase the status line.
		case '\r':
		    call printf ("\n")

		# Print help.
		case '?':
		    ;

		# Process the next image.
		case 'n':
		    if (im == NULL) {
			call xp_stats (xp, IMAGE, Memc[imname], SZ_FNAME)
			if (Memc[imname] == EOS)
		            call printf (
			    "Warning: The image is undefined\n")
			else
		            call printf (
			    "Warning: Cannot open image (%s)\n")
			        call pargstr (Memc[imname])
		    } else if (imno < imtlen (imlist)) {
			NEWIMAGE(pstatus) = YES
			if (olno < fntlenb (objlist))
			    NEWLIST(pstatus) = YES
			if (rlno < fntlenb (reslist))
			    NEWRESULTS(pstatus) = YES
			if (glno < fntlenb (greslist))
			    NEWRESULTS(pstatus) = YES
		        break
		    } else {
		        call printf ("Warning: The image list is at EOF\n")
		    }

		# Process the previous image.
		case 'p':
		    if (im == NULL) {
			call xp_stats (xp, IMAGE, Memc[imname], SZ_FNAME)
			if (Memc[imname] == EOS)
		            call printf (
			    "Warning: The image is undefined\n")
			else
		            call printf (
			    "Warning: Cannot open image (%s)\n")
			        call pargstr (Memc[imname])
		    } else if (imno > 1) {
		        imno = imno - 2
			NEWIMAGE(pstatus) = YES
			if (olno > 1) {
			    olno = olno - 2
			    NEWLIST(pstatus) = YES
			}
			if (rlno > 1) {
			    rlno = rlno - 2
			    NEWRESULTS(pstatus) = YES
			}
			if (glno > 1) {
			    glno = glno - 2
			    NEWRESULTS(pstatus) = YES
			}
		        break
		    } else {
		        call printf ("Warning: The image list is at BOF\n")
		    }

		# Redisplay the image.
		case 'i':
		    if (im == NULL) {
			call xp_stats (xp, IMAGE, Memc[imname], SZ_FNAME)
			if (Memc[imname] == EOS)
		            call printf (
			    "Warning: The image is undefined\n")
			else
		            call printf (
			    "Warning: Cannot open image (%s)\n")
			        call pargstr (Memc[imname])
		    } else {
	    	        call xp_display (gd, xp, im, 1, IM_LEN(im,1), 1,
		            IM_LEN(im,2), IMAGE_DISPLAY_WCS, IMAGE_DISPLAY_WCS)
			call xp_keyset (im, xp)
		    }

		# Display the header of the image.
		case 'h':
		    if (im != NULL) {
	    	        call xp_pimheader (gd, im)
		    } else {
			call xp_stats (xp, IMAGE, Memc[imname], SZ_FNAME)
			if (Memc[imname] == EOS)
		            call printf (
			    "Warning: The image is undefined\n")
			else
		            call printf (
			    "Warning: Cannot open image (%s)\n")
			        call pargstr (Memc[imname])
		    }

		# Process the next coordinate list.
		case ']':
		    if (olno < fntlenb (objlist)) {
			NEWLIST(pstatus) = YES
		        break
		    }

		# Process the previous coordinate list.
		case '[':
		    if (olno > 1) {
			olno = olno - 2
			NEWLIST(pstatus) = YES
			break
		    }

		# Create, examine, and save the objects symbol tables.
		case '^', 'f', 'b', '~', 'm', 'e', '@', 'a', 'd', 'u', 'z', 'r',
		    'l', 'w':
		    symbol = xp_dmeas (gd, xp, im, ol, rl, key, wx, wy, NO)

		# Do a quick model fit and print results.
		case 'x':
		    iupdate = xp_pmodfit (gd, xp, im, wx, wy, rl, 4, NO, NO)

		# Do a quick model fit and plot results.
		case 'y':
		    iupdate = xp_pmodfit (gd, xp, im, wx, wy, rl, 4, YES, NO)

		# Do a quick model fit, plot and interact with the plot.
		case 'Y':
		    iupdate = xp_pmodfit (gd, xp, im, wx, wy, rl, 4, YES, YES)

		# Make an objects list.
		case ' ':
		    if (rl != NULL && LOGRESULTS(pstatus) == YES) {
			nobjs = xp_afind (gd, im, rl, xp, NO)
			call printf (
			    "Saved %d objects in results file %s\n")
			    call pargi (nobjs)
			    call pargstr (Memc[rlname])
		    } else {
			nobjs = xp_afind (gd, im, NULL, xp, NO)
			call printf (
			    "Detected %d objects in image %s\n")
			    call pargi (nobjs)
			    call pargstr (Memc[imname])
		    }
		    symbol = NULL

                # Display a zoomed-up subraster around the measured star.
                case 'I':
                    if (im != NULL) {
                        if (symbol == NULL) {
                            call xp_cpdisplay (gd, xp, im, wx, wy, INDEFR,
			        IMAGE_DISPLAY_WCS, IMAGE_DISPLAY_WCS, YES, NO)
                        } else {
                            call xp_ocpdisplay (gd, xp, im, XP_OXINIT(symbol),
			        XP_OYINIT(symbol), INDEFR, IMAGE_DISPLAY_WCS,
				IMAGE_DISPLAY_WCS, YES, NO)
                        }
                    } else {
                        call xp_stats (xp, IMAGE, Memc[imname], SZ_FNAME)
                        call printf ("Warning: Cannot access image (%s)\n")
                            call pargstr (Memc[imname])
                    }

                # Draw a contour plot of a subraster around the measured star.
                case 'c':
                    if (im != NULL) {
                        if (symbol == NULL)
                            call xp_cpplot (gd, xp, im, wx, wy, INDEFR,
			        IMAGE_DISPLAY_WCS, IMAGE_DISPLAY_WCS)
                        else
                            call xp_ocpplot (gd, xp, im, XP_OXINIT(symbol),
                                XP_OYINIT(symbol), INDEFR,
				IMAGE_DISPLAY_WCS, IMAGE_DISPLAY_WCS)
                    } else {
                        call xp_stats (xp, IMAGE, Memc[imname], SZ_FNAME)
                        call printf ("Warning: Cannot access image (%s)\n")
                            call pargstr (Memc[imname])
                    }

                # Draw a contour plot of a subraster around the measured star.
                case 's':
                    if (im != NULL) {
                        if (symbol == NULL)
                            call xp_asplot (gd, xp, im, wx, wy, INDEFR,
			        IMAGE_DISPLAY_WCS, IMAGE_DISPLAY_WCS)
                        else
                            call xp_oasplot (gd, xp, im, XP_OXINIT(symbol),
                                XP_OYINIT(symbol), INDEFR, IMAGE_DISPLAY_WCS,
				IMAGE_DISPLAY_WCS)
                    } else {
                        call xp_stats (xp, IMAGE, Memc[imname], SZ_FNAME)
                        call printf ("Warning: Cannot access image (%s)\n")
                            call pargstr (Memc[imname])
                    }

                # Overlay a contour plot on a display of subraster around
                # the measured star.
                case 'O':
                    if (im != NULL) {
                        if (symbol == NULL) {
                            call xp_cpdisplay (gd, xp, im, wx, wy, INDEFR,
			        IMAGE_DISPLAY_WCS, IMAGE_DISPLAY_WCS, YES, YES)
                        } else {
                            call xp_ocpdisplay (gd, xp, im, XP_OXINIT(symbol),
			        XP_OYINIT(symbol), INDEFR, IMAGE_DISPLAY_WCS,
				IMAGE_DISPLAY_WCS, YES, YES)
                        }
                    } else {
                        call xp_stats (xp, IMAGE, Memc[imname], SZ_FNAME)
                        call printf ("Warning: Cannot access image (%s)\n")
                            call pargstr (Memc[imname])
                    }

		# Process a colon command.
		case ':':
		    symbol = xp_xdcolon (gd, xp, dirlist, imlist, im,
		        objlist, ol, reslist, rl, greslist, gl, Memc[cmd],
			symbol)

		# Unknown or ambiguous command.
		default:
		    call printf ("Ambiguous or undefined keystroke command\n")
		}

		# Display the new image.
		if (NEWIMAGE(pstatus) == YES || NEWLIST(pstatus) == YES ||
		    NEWRESULTS(pstatus) == YES) {
		    if (NEWIMAGE(pstatus) == YES) {
		        call xp_stats (xp, IMAGE, Memc[imname], SZ_FNAME)
		        imno = xp_stati (xp, IMNUMBER)
		        if (im != NULL) {
		            call xp_display (gd, xp, im, 1, IM_LEN(im,1), 1,
		                IM_LEN(im,2), IMAGE_DISPLAY_WCS,
				IMAGE_DISPLAY_WCS)
		        } else {
			    call gclear (gd)
			    call gflush (gd)
		            call printf (
			    "Warning: Cannot open image (%s)\n")
			        call pargstr (Memc[imname])
		        }
			call xp_keyset (im, xp)
		    }
		    if (NEWLIST(pstatus) == YES) {
			OBJNO(pstatus) = 0
		        call xp_stats (xp, OBJECTS, Memc[olname], SZ_FNAME)
		        olno = xp_stati (xp, OFNUMBER)
		        if (ol != NULL) {
			    ;
		        } else if (fntlenb (objlist) > 0) {
		            call printf (
			    "Warning: Cannot open object list (%s)\n")
			        call pargstr (Memc[olname])
		        }
		    }
		    if (NEWRESULTS(pstatus) == YES) {
		        call xp_stats (xp, RESULTS, Memc[rlname], SZ_FNAME)
		        rlno = xp_stati (xp, RFNUMBER)
		        if (rl != NULL) {
                            if (SEQNO(pstatus) == 0)
				call xp_whfind (xp, rl, "xdisplay")
                            call xp_whiminfo (xp, rl)
                            if (SEQNO(pstatus) == 0)
				call xp_xfbnr (xp, rl)
		        } else if (fntlenb (reslist) > 0) {
		            call printf (
			    "Warning: Cannot open results list (%s)\n")
			        call pargstr (Memc[rlname])
			}
                        call xp_stats (xp, GRESULTS, Memc[glname], SZ_FNAME)
                        glno = xp_stati (xp, GFNUMBER)
                        if (gl != NULL) {
                            ;
                        } else if (fntlenb (greslist) > 0) {
                            call printf (
                            "Warning: Cannot open results list (%s)\n")
                                call pargstr (Memc[glname])
                        }
		    }
		    if (NEWIMAGE(pstatus) == YES) {
		        owx = INDEFR; owy = INDEFR
		    } else if (NEWLIST(pstatus) == YES) {
		        owx = wx; owy = wy
		    }
		    NEWIMAGE(pstatus) = NO
		    NEWLIST(pstatus) = NO
		    NEWRESULTS(pstatus) = NO
		} else {
		    owx = wx; owy = wy
		}
	    }

noninteractive_
	    if (gd == NULL) {
		nobjs = xp_afind (NULL, im, rl, xp, verbose)
	    }

	    # Increment the results file list counter.
	    if (imno == EOF || NEWRESULTS(pstatus) == YES) {
	        if (rl != NULL) {
		    call flush (rl)
		    if (fstati (rl, F_FILESIZE) == 0 || SEQNO(pstatus) == 0) {
			call fstats (rl, F_FILENAME, Memc[str], SZ_FNAME)
			call close (rl)
			call delete (Memc[str])
		    } else {
		        call close (rl)
                        lsymbol = stfind (xp_statp(xp, SEQNOLIST), Memc[rlname])
                        if (lsymbol != NULL)
                            XP_MAXSEQNO(lsymbol) = SEQNO(pstatus)
                        else {
                            lsymbol = stenter (xp_statp(xp, SEQNOLIST),
                                Memc[rlname], LEN_SEQNOLIST_STRUCT)
                            XP_MAXSEQNO(lsymbol) = SEQNO(pstatus)
                        }
		    }
		    rl = NULL
	        }
	        if ((imno != EOF) && (rlno < fntlenb (reslist))) 
		    rlno = rlno + 1
                if (gl != NULL) {
                    call flush (gl)
                    if (fstati (gl, F_FILESIZE) == 0 || SEQNO(pstatus) == 0) {
                        call fstats (gl, F_FILENAME, Memc[str], SZ_FNAME)
                        call close (gl)
                        call delete (Memc[str])
                    } else
                        call close (gl)
                    gl = NULL
                }
                if ((imno != EOF) && (glno < fntlenb (greslist)))
                    glno = glno + 1
	    }

	    # Increment the object list counter.
	    if (imno == EOF || NEWLIST(pstatus) == YES) {
	        if (ol != NULL) {
		    call close (ol)
		    ol = NULL
	        }
	        if ((imno != EOF) && (olno < fntlenb (objlist))) 
		    olno = olno + 1
	    }

	    # Increment the image list counter.
	    if (imno == EOF || NEWIMAGE(pstatus) == YES) {
	        if (im != NULL) {
	            call imunmap (im)
		    im = NULL
	        }
	        if ((imno != EOF) && (imno < imtlen (imlist)))
	            imno = imno + 1
		else if (gd == NULL)
		    imno = EOF
	    }

	} until (imno == EOF)

	# Update the algorithm parameters.
	if (update == YES)
	    call xp_pxdpars (xp)

	# Reset the current directory
	call xp_stats (xp, STARTDIR, Memc[cmd], SZ_LINE)
	call fchdir (Memc[cmd])

	call xp_xdfree (xp)
	call fntclsb (greslist)
	call fntclsb (reslist)
	call fntclsb (objlist)
	call imtclose (imlist)
	call fntclsb (dirlist)
	if (gd != NULL)
	    call gclose (gd)

	call sfree (sp)
end
