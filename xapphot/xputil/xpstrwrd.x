# XP_STRWRD -- Search a dictionary string for a given string index number.
# This is the opposite function of strdic(), that returns the index for
# given string.  The entries in the dictionary string are separated by
# a delimiter character which is the first character of the dictionary
# string.  The index of the string found is returned as the function value.
# Otherwise, if there is no string for that index, a zero is returned.

int procedure xp_strwrd (index, outstr, maxch, dict)

int	index			#I the string index
char	outstr[ARB]		#O the output string as found in dictionary
int	maxch			#I the maximum length of output string
char	dict[ARB]		#I the dictionary string

int	i, len, start, count

int	strlen()

begin
	# Clear output string
	outstr[1] = EOS

	# Return if the dictionary is not long enough
	if (dict[1] == EOS)
	    return (0)

	# Initialize counters
	count = 1
	len   = strlen (dict)

	# Search the dictionary string. This loop only terminates
	# successfully if the index is found. Otherwise the procedure
	# returns with and error condition.
	for (start = 2; count < index; start = start + 1) {
	    if (dict[start] == dict[1])
		count = count + 1
	    if (start == len)
		return (0)
	}

	# Extract the output string from the dictionary
	for (i = start; dict[i] != EOS && dict[i] != dict[1]; i = i + 1) {
	    if (i - start + 1 > maxch)
		break
	    outstr[i - start + 1] = dict[i]
	}
	outstr[i - start + 1] = EOS

	# Return index for output string
	return (count)
end