include "pstools.h"


# PS_CLOSE -- Close the struct and free memory.

procedure ps_close (ps)

pointer	ps					#i PSTOOLS descriptor

begin
	# Write the page trailer.
	call ps_trailer(ps)

	call mfree (PS_HLE(ps), TY_CHAR)
	call mfree (PS_HCE(ps), TY_CHAR)
	call mfree (PS_HRE(ps), TY_CHAR)
	call mfree (PS_FLE(ps), TY_CHAR)
	call mfree (PS_FCE(ps), TY_CHAR)
	call mfree (PS_FRE(ps), TY_CHAR)
	call mfree (PS_WBPTR(ps), TY_CHAR)

	call mfree (ps, TY_STRUCT)
end
