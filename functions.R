###
# Long functions for the analysis
# Zane
# 2025-05-21
# Most of these functions are copied from another source which will be cited.
###

# Functions from agdist repo sequence-alignment.R code ####
# https://github.com/ahgroup/influenza-antigenic-distance-data/blob/33d9b7607176207cc37de23df1c7bfa84705673b/R/sequence-alignment.R
clean_sequence_data <- function(raw_sequence_data) {
	out <-
		raw_sequence_data |>
		# clean up the column names
		janitor::clean_names() |>
		# Combine the accession number columns
		tidyr::unite(
			col = "accession_no",
			c(protein_sequence_source, protein_sequence_accession_number),
			sep = ": "
		) |>
		# Cleanup
		dplyr::mutate(
			# Convert strain name to short format -- needed to do this once hgp
			# dependence was removed from this code, all downstream code depends
			# on the short name beinged named "strain_name" and it's not worth the
			# effort to change that.
			full_name = strain_name,
			strain_name = short_name,
			# Remove any newlines from the sequences
			protein_sequence = stringr::str_remove_all(protein_sequence, "\\s"),
			# Ensure all sequences are lowercase
			protein_sequence = stringr::str_to_lower(protein_sequence),
			# Remove any excess whitespace from all columns
			dplyr::across(
				dplyr::where(is.character),
				stringr::str_squish
			),
			# Add a column with the sequence length
			sequence_length = nchar(protein_sequence)
		) |>
		# Remove the duplicate name column
		dplyr::select(-short_name)
	
	return(out)
}

align_sequence_data <- function(processed_sequence_data) {
	# we'll also give an easier to access name to the strain names data
	unaligned_sequences <-
		processed_sequence_data |>
		# Remove the accession_no column that we don't need for analysis
		dplyr::select(-accession_no) |>
		# Add a subtype column that we can use later, and also add a column for
		# the longer strain name.
		dplyr::mutate(
			analysis_name = hgp::replace_strain_names(
				strain_name, from = "short", to = "analysis"
			),
			subtype = hgp::replace_strain_names(
				strain_name, from = "short", to = "subtype"
			),
			.after = strain_name
		) |>
		# Remove the overall strain, we don't need that here
		dplyr::filter(subtype != "") |>
		dplyr::mutate(dplyr::across(dplyr::where(is.factor), forcats::fct_drop)) |>
		# Make a new type/subtype variable that allows us to separate the A(H1N1)
		# and A(H3N2) strains while keeping all B strains together for alignment.
		dplyr::mutate(
			type_subtype = dplyr::case_when(
				strain_type == "B" ~ "B",
				strain_subtype == "H1N1" ~ "A(H1N1)",
				strain_subtype == "H3N2" ~ "A(H3N2)"
			)
		)
	
	# Remove signal peptides -- see these sources for signal peptide length:
	# TODO add these
	# TODO we currently don't remove the signal peptides, do we even actually
	# need to do this?
	signal_peptide_end <- dplyr::case_match(
		unaligned_sequences$type_subtype,
		"A(H1N1)" ~ 18L,
		"A(H3N2)" ~ 17L,
		"B" ~ 15L,
		.default = NA_integer_
	)
	
	# Next we'll manually align the sequences. We tried using a multiple alignment
	# algorithm as a shortcut, but they didn't correctly align all of the
	# conserved indel positions. So we used that to align the starting locations
	# of the partial length sequences, and then we used literature knowledge to
	# manually align the indels in the sequences.
	seqs_aligned <-
		unaligned_sequences |>
		dplyr::mutate(
			# Align the beginnings of the partial length sequences
			# This is based on a temporary MSA using MUSCLE -- we aligned the
			# H1N1, H3N2, and B strains separately, and then extracted the starting
			# positions of the partial length sequences from the MSA.
			protein_sequence = dplyr::case_when(
				strain_name == "MI/85" ~ paste0("m", strrep("x", 15), protein_sequence),
				strain_name == "Sing/64" ~ paste0("m", strrep("x", 3), protein_sequence),
				strain_name == "Sich/99" ~ paste0("m", strrep("x", 7), protein_sequence),
				TRUE ~ protein_sequence
			),
			# Recompute the sequence length
			intermediate_length = nchar(protein_sequence),
			# Align the gaps in the H1N1 sequences -- all of the H1N1 sequences are of
			# either length 566 (CA/09-like and SC/18-like lineages) or length 565
			# (non-pdm-like lineages). The sequences of length 565 need to have a
			# gap inserted at position 130+18 (i.e. the 130th residue after the
			# 18 residue signal peptide). Need to find a citation for this since I
			# just got it from Amanda. This regular expression puts a - after the
			# first 129 + 18 characters so that the dash ends up at spot 130 + 18.
			protein_sequence = dplyr::if_else(
				sequence_length == 565 & subtype == "H1N1",
				# This regular expression adds a dash at spot 148 including SP.
				gsub('^(.{147})(.*)$', '\\1-\\2', protein_sequence),
				protein_sequence
			),
			# Now we need to align the indels in the flu B sequences, so that all of
			# the flu B sequences are of length 585 residues. We need them to be
			# all aligned so we can correctly calculate the antigenic distance across
			# the flu B lineages.
			# Gaps are included based on Lee40 numbering following this paper:
			# https://pmc.ncbi.nlm.nih.gov/articles/PMC104260/
			# https://pmc.ncbi.nlm.nih.gov/articles/PMC7565853/
			protein_sequence = dplyr::case_when(
				# 1DEL lineages (all pre-split and yamagata, 2 victoria)
				# Maybe more accurate that the victoria lineages have one insertion
				strain_type == "B" & intermediate_length == 584 ~
					gsub('^(.{162})(.*)$', '\\1-\\2', protein_sequence),
				# 2DEL lineages of B-victoria and 1 yamagata
				strain_type == "B" & intermediate_length == 583 ~
					gsub('^(.{162})(.*)$', '\\1-\\2',
							 gsub('^(.{164})(.*)$', '\\1-\\2', protein_sequence)),
				# 3DEL lineages of B-victoria
				strain_type == "B" & intermediate_length == 582 ~
					gsub('^(.{162})(.*)$', '\\1-\\2',
							 gsub('^(.{163})(.*)$', '\\1-\\2',
							 		 gsub('^(.{164})(.*)$', '\\1-\\2', protein_sequence))),
				TRUE ~ protein_sequence
			),
			# Next, add "X" characters to the end of the partial length sequences.
			protein_sequence = dplyr::case_when(
				# H1 and H3 should have final length 566
				strain_type == "A" ~
					stringr::str_pad(
						protein_sequence,
						width = 566,
						side = "right",
						pad = "x"
					),
				# B should have final length 585
				strain_type == "B" ~
					stringr::str_pad(
						protein_sequence,
						width = 585,
						side = "right",
						pad = "x"
					),
				TRUE ~ protein_sequence
			),
			# Finally compute the length of the aligned sequence so we can make sure
			# that we've finished all of the alignment steps.
			aligned_sequence_length = nchar(protein_sequence)
		) |>
		# Rename the plain sequence_length variable
		dplyr::rename(unaligned_sequence_length = sequence_length) |>
		# And get rid of the intermediate one
		dplyr::select(-intermediate_length)
	
	return(seqs_aligned)
}

# From influenza antigenic distance utils.R file
# https://github.com/ahgroup/influenza-antigenic-distance-data/blob/33d9b7607176207cc37de23df1c7bfa84705673b/R/utils.R

# Function for excluding ambiguous residues coded as X in two sequences
# This will do listwise removal of the specified ambiguous residues across
# all elements of the vector or list 'seqs'.
# If you want to do pairwise removal, you should invoke this on only two
# different sequences.
remove_ambiguous_residues <- function(seqs, ambiguous_residues = c("xX?")) {
	
	# First validate the seqs arugment to make sure it's a character vector
	# or a list of character vectors
	validate_sequence_input_form(seqs)
	
	# Validate the ambiguous residues argument and split it into a vector
	if (is.null(ambiguous_residues) || is.na(ambiguous_residues)) {
		stop("'ambiguous_residues' argument is NULL or NA.")
	} else if (!is.character(ambiguous_residues) || length(ambiguous_residues) > 1L) {
		stop("'ambiguous_residues' argument should be a length 1 character vector.")
	} else if (nchar(ambiguous_residues) == 0) {
		# If the ambiguous residues argument is the empty string "", return the
		# original seqs argument (after validation). This provides an easy way
		# to skip removal if needed.
		return(seqs)
	} else {
		ambiguous_residues <- strsplit(ambiguous_residues, "")[[1]]
	}
	
	# Now if seqs is a character vector it needs to be split into a list of
	# character vectors where all elements are length one.
	if (is.character(seqs)) {
		seq_split <- strsplit(seqs, "")
	} else {
		seq_split <- seqs
	}
	
	# Find out which elements of seqs are ambiguous.
	ambiguous_indices <- sapply(
		seq_split,
		\(x) which(x %in% ambiguous_residues)
	)
	# Now combine them into one list of all the ambiguous sites.
	all_ambiguous_indices <-
		Reduce(union, ambiguous_indices, init = integer(0)) |>
		sort()
	
	# Remove the ambiguous ones from each element of seq_split -- if the vector
	# of indices is empty though we have to handle this separately.
	if (length(all_ambiguous_indices) > 0) {
		seq_unambiguous <-
			lapply(seq_split, \(x) x[-all_ambiguous_indices])
	} else {
		seq_unambiguous <- seq_split
	}
	
	# If the input vector was a character, return a character vector -- we assume
	# the user wants output in the same format as the input.
	if (is.character(seqs)) {
		seqs_out <- sapply(seq_unambiguous, \(x) paste0(x, collapse = ""))
	} else {
		seqs_out <- seq_unambiguous
	}
	
	return(seqs_out)
}

# Function to run a stringdist distance function while removing ambiguous
# residues
dist_string <- function(seqs, ambiguous_residues = c("xX?"), ...) {
	# Set up an empty matrix to hold results
	res <- matrix(
		nrow = length(seqs),
		ncol = length(seqs),
		dimnames = list(names(seqs), names(seqs))
	)
	
	# Calculate the lower triangle of the distance matrix for all unique
	# combinations of two strains
	# This invokes stringdist::stringdist and the ... arguments are only
	# evaluated in this context.
	for (i in 2:nrow(res)) {
		for (j in 1:(i-1)) {
			seqs_clean <- remove_ambiguous_residues(
				c(seqs[[i]], seqs[[j]]),
				ambiguous_residues
			)
			res[[i, j]] <- stringdist::stringdist(
				a = seqs_clean[[1]], b = seqs_clean[[2]],
				...
			) / nchar(seqs_clean[[1]])
		}
	}
	
	# Set the diagonal to zero
	diag(res) <- 0
	
	# Quickly fill in the upper triangle
	out <- res |> as.dist() |> as.matrix()
	
	return(out)
}

# This function checks if a sequence argument is either a character vector or
# a list of pre-splitted sequences, i.e. a list of character vectors that
# are all length one.
validate_sequence_input_form <- function(seqs, require_alignment = TRUE) {
	if(is.character(seqs)) {
		unique_lengths <- unique(nchar(seqs))
		if(length(unique_lengths) != 1) {
			stop("seqs is a character vector, but sequences are not aligned.")
		} else {
			return(TRUE)
		}
	}
	
	if(!is.list(seqs)) {
		stop("'seqs' argument is not a character vector or a list.")
	}
	
	if(!all(sapply(seqs, is.character))) {
		stop("'seqs' was a list, but at least one element was not a character vector.")
	}
	
	if((length(unique(sapply(seqs, length))) != 1) && isTRUE(require_alignment)) {
		stop(paste0(
			"'seqs' was a list of character vectors, but they are not aligned.",
			" They should all be the same length."
		))
	}
	
	if(!all(sapply(seqs, \(x) all(lapply(x, length) == 1)))) {
		stop(paste0(
			"'seqs' was a list of character vectors, but at least one element was ",
			"longer than length one."
		))
	}
	
	return(TRUE)
}

# Utility function for subsetting characters in a single string (in R language,
# getting a substring by position for a character vector of length 1)
extract_string_chars <- function(str, pos) {
	sub <- charToRaw(str)[pos]
	return(rawToChar(sub))
}

# Utility function to set the row and column names of a square matrix to be
# the same, in one function call.
set_square_matrix_dimnames <- function(matrix, names) {
	names <- as.character(names)
	out <- matrix
	colnames(out) <- names
	rownames(out) <- names
	
	return(out)
}

# Convert a distance matrix to tidy data format
tidy_dist_mat <- function(d) {
	out <- d |>
		tibble::as_tibble(rownames = "Strain1") |>
		tidyr::pivot_longer(
			cols = -Strain1,
			names_to = "Strain2",
			values_to = "d"
		) |>
		# Order variable factors
		dplyr::mutate(
			Strain1 = forcats::fct_inorder(as.character(Strain1)),
			Strain2 = forcats::fct_inorder(as.character(Strain2)) |> forcats::fct_rev()
		)
	
	return(out)
}

# From peptiope calculation script
# https://github.com/ahgroup/influenza-antigenic-distance-data/blob/33d9b7607176207cc37de23df1c7bfa84705673b/R/pepitope.R

# First we need to create a reference data frame of p-epitope site locations
get_pepitope_sites <- function(
		subtype,
		sites = c('a', 'b', 'c', 'd', 'e', 'all_epitopes'),
		harmonize_b_lineages = TRUE
) {
	
	# standardize subtype
	subtype_length <- nchar(subtype)
	subtype_san <-
		# If the subtype has parentheses, assume it is in the format A(HXNY) or
		# B(lineage), otherwise assume it is in the format HXNY or lineage.
		# If there's parantheses, first we remove them.
		dplyr::if_else(
			stringr::str_detect(subtype, "\\(", negate = TRUE),
			subtype,
			stringr::str_extract(subtype, "\\((.{4})\\)", group = 1)
		) |>
		# Then we sanitize the input so we don't have to care about capitalization
		# and potentially fix misspellings.
		stringr::str_to_lower() |>
		substr(1, min(2, subtype_length))
	
	# Validate the harmonize_b_lineages argument
	if (!isTRUE(harmonize_b_lineages) && !isFALSE(harmonize_b_lineages)) {
		stop(paste0(
			"the 'harmonize_b_lineages' argument should be TRUE or FALSE, not ",
			harmonize_b_lineages
		))
	}
	
	# Validate the site argument
	if (!all(sites %in% c('a', 'b', 'c', 'd', 'e', 'all_epitopes'))) {
		stop(paste0(
			"sites should be a character vector containing one or more of the ",
			"following allowed entries: 'a', 'b', 'c', 'd', 'e', or 'all_epitopes', ",
			"instead of : ",
			sites
		))
	}
	
	# Check if the B subtypes should be combined
	if (
		(isTRUE(harmonize_b_lineages) && (subtype_san %in% c("vi", "ya", "pr"))) ||
		(startsWith(subtype_san, "b"))
	) {
		subtype_san <- "b_unified"
	}
	
	# The sites are different for H1 and H3. These are based on Amanda's
	# code where she referenced the actual papers.
	if (subtype_san == "h1") {
		#full <- 1:326
		
		site_a <- c(
			118, 120:122, 126:129, 132:135, 137, 139:143, 146, 147, 149, 165, 252,
			253
		)
		
		site_b <- c(
			124, 125, 152:157, 160, 162, 183:187, 189:191, 193:196
		)
		
		site_c <- c(
			34:38, 40, 41, 43:45, 269:274, 276:278, 283, 288, 292, 295, 297, 298,
			302, 303, 305:310
		)
		
		site_d <- c(
			89, 94:96, 113, 117, 163, 164, 166:174, 176, 179, 198, 200, 202, 204:216,
			222:227, 235, 237, 241, 243:245
		)
		
		site_e <- c(
			47, 48, 50, 51, 53, 54, 56:58, 66, 68:75, 78:80, 82:86, 102, 257:261,
			263, 267
		)
		
	} else if (subtype_san == "h3") {
		#full <- 1:328
		
		site_a <- c(
			122, 124, 126, 130:133, 135, 137:138, 140, 142:146, 150, 152,
			168
		)
		
		site_b <- c(
			128, 129, 155:160, 163, 165, 186:190, 192:194, 196:198
		)
		
		site_c <- c(
			44:48, 50, 51, 53, 54, 273, 275:276, 278:280, 294, 297, 299, 300,
			304:305, 307:312
		)
		
		site_d <- c(
			96, 102, 103, 117, 121, 167, 170:177, 179, 182, 201, 203, 207:209,
			212:219, 226:230, 238, 240, 242, 244, 246:248
		)
		
		site_e <- c(
			57, 59, 62, 63, 67, 75, 78, 80:83, 86:88, 91, 92, 94, 109, 260:262,
			265
		)
	} else if (subtype_san == "b_unified") {
		# Case for combining the two b lineages together
		site_a <- c(
			121, 122, 123, 125, 126, 134, 135, 136, 137, 139, 141, 142, 144, 146,
			147, 148, 149, 150, 151, 155, 157, 176, 177
		)
		
		site_b <- c(
			127, 129, 133, 160, 161, 162, 163, 164, 165, 166, 167, 168, 171, 172,
			173, 174, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206,
			207, 208, 209
		)
		
		site_c <- c(
			34, 35, 36, 37, 38, 39, 40, 288, 289, 290, 291, 292, 293, 294, 308, 309,
			314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327
		)
		
		site_d <- c(
			93, 101, 102, 116, 120, 175, 176, 178, 179, 180, 181, 182, 183, 184, 185,
			186, 187, 188, 189, 190, 211, 212, 213, 214, 217, 218, 219, 220, 222,
			223, 224, 225, 226, 227, 228, 229, 230, 232, 233, 241, 242, 243, 244,
			245, 246, 253, 254, 255, 256, 257, 258
		)
		
		site_e <- c(
			42, 44, 48, 56, 58, 59, 63, 71, 73, 75, 77, 78, 79, 80, 83, 84, 85, 88,
			89, 91, 108, 272, 273, 275, 276, 277, 279, 280
		)
		
	} else if (subtype_san == "vi") {
		# Case for victoria lineage with B subtypes lineage
		site_a <- c(
			121, 122, 123, 125, 126, 134, 135, 136, 137, 139, 141, 142, 144, 146,
			147, 148, 149, 150, 151, 155, 157, 177
		)
		
		site_b <- c(
			127, 129, 133, 160, 161, 162, 163, 164, 165, 166, 168, 172, 174,
			196, 197, 198, 199, 200, 202, 203, 204, 206, 207, 208, 209
		)
		
		site_c <- c(
			34, 35, 36, 37, 38, 39, 40, 289, 291, 292, 293, 294, 309, 315,
			317, 318, 320, 321, 323, 324, 325, 326, 327
		)
		
		site_d <- c(
			93, 101, 102, 116, 120, 176, 179, 180, 182, 183, 184, 185, 186,
			187, 188, 190, 212, 214, 218, 219, 220, 223, 224, 225, 226,
			227, 228, 229, 230, 233, 242, 243, 244, 245, 246, 254, 255,
			256, 257, 258
		)
		
		site_e <- c(
			42, 44, 48, 56, 58, 59, 63, 71, 73, 75, 77, 78, 79, 80, 83, 84, 85,
			88, 89, 91, 108, 273, 276, 277, 280
		)
		
	} else if (subtype_san %in% c("ya", "pr")) {
		# Case for yamagata and pre-divergence lineages with B lineages separate
		site_a <- c(
			121, 122, 123, 125, 126, 134, 135, 136, 137, 139, 141, 142, 144,
			146, 147, 148, 149, 150, 151, 155, 157, 176
		)
		
		site_b <- c(
			127, 129, 133, 160, 161, 162, 163, 164, 165, 167, 171, 173, 195,
			196, 197, 198, 199, 201, 202, 203, 205, 206, 207, 208
		)
		
		site_c <- c(
			34, 35, 36, 37, 38, 39, 40, 288, 290, 291, 292, 293, 308, 314,
			316, 317, 319, 320, 322, 323, 324, 325, 326
		)
		
		site_d <- c(
			93, 101, 102, 116, 120, 175, 178, 179, 181, 182, 183, 184, 185,
			186, 187, 189, 211, 213, 217, 218, 219, 222, 223, 224,
			225, 226, 227, 228, 229, 232, 241, 242, 243, 244, 245, 253, 254,
			255, 256, 257
		)
		
		site_e <- c(
			42, 44, 48, 56, 58, 59, 63, 71, 73, 75, 77, 78, 79, 80, 83, 84, 85,
			88, 89, 91, 108, 272, 275, 276, 279
		)
	} else {
		# If the subtype didn't match one of the above arguments, then it failed
		# validation.
		stop(paste0(
			"subtype should be one of the following: h1n1, h3n2, victoria, ",
			"yamagata, or presplit. Capitalization doesn't matter and abbreviations ",
			"are allowed. Instead you specified: ",
			subtype_san,
			" which is not supported."
		))
	}
	
	# Make a vector containing all epitope residues
	all <- c(site_a, site_b, site_c, site_d, site_e) |>
		unique() |>
		sort()
	
	# Make a list containing all the vectors
	res <- list(site_a, site_b, site_c, site_d, site_e, all)
	names(res) <- c('a', 'b', 'c', 'd', 'e', 'all_epitopes')
	
	# Now return the sites that were requested before
	out <- res[sites]
	return(out)
}

# Compute the p-epitope distance between two strings
pepitope <- function(seq_1, seq_2, subtype, mode = "dominant", harmonize_b_lineages = TRUE, ambiguous_residues) {
	# Get the numbers for the residues in each epitope
	p_epi_sites <- get_pepitope_sites(
		subtype,
		sites = c('a', 'b', 'c', 'd', 'e'),
		harmonize_b_lineages = TRUE
	)
	
	epi_dists <-
		purrr::map_dbl(
			p_epi_sites,
			# This function looks kind of weird because doing this calculation is
			# really annoying to make sure the Hamming dist is the same as the
			# overall Hamming dist.
			# First we subset the two strings to only contain residues from the
			# current epitope.
			\(current_sites) list(
				"a" = extract_string_chars(seq_1, current_sites),
				"b" = extract_string_chars(seq_2, current_sites),
				method = "hamming"
			) |>
				do.call(what = stringdist::stringdist)
		)
	
	# Normalize the distances by the length of each site
	epi_dists_raw <- epi_dists
	epi_dists <- epi_dists / lengths(p_epi_sites)
	
	# Now decide whether to return all epitopes or not
	if (mode %in% c("dominant", "max", "maximum")) {
		# Dominant p-epitope = maximum p-epitope
		# See Gupta 2006 (PMID 16460844).
		p_epi <- max(epi_dists)
	} else if (mode == "anderson") {
		# In Anderson 2018 (PMID 29433425) they take the average of all of the
		# epitope differences. They call this "p-all-epitope" and it is equivalent
		# in practice but provided here just in case. We leave out the scaling by
		# 20 but otherwise the formula is reproduced exactly from the paper.
		p_epi <- mean(epi_dists / lengths(p_epi_sites))
	} else if (mode  %in% c("all", "average", "mean")) {
		# THis is the original definition of p-all-epitope from Pan et al (2010),
		# PMID: 21123189. They basically look at the Hamming distance but only
		# considering the sites that are in the epitopes.
		# This is the same as the mean of the epitope distances, but we use
		# the exact formula of the paper for posterity.
		p_epi <- sum(epi_dists) / sum(lengths(p_epi_sites))
	} else if (mode == "median") {
		# Median instead of mean because Andreas likes the median.
		p_epi <- median(epi_dists)
	} else if (is.null(mode) | is.na(mode) | mode == "") {
		# If there's no summary mode specified, return the entire vector of
		# differences.
		p_epi <- epi_dists
	} else {
		rlang::warn("p-Epitope mode specified with an option that doesn't exist.")
		p_epi <- epi_dists
	}
	return(p_epi)
}

# Compute a p-epitope distance matrix for all sequences in a character vector
dist_pepi <- function(seqs, subtype, mode = "dominant", harmonize_b_lineages = TRUE, ambiguous_residues = c("xX?")) {
	# Set up an empty matrix to hold results
	res <- matrix(
		nrow = length(seqs),
		ncol = length(seqs),
		dimnames = list(names(seqs), names(seqs))
	)
	
	# Calculate the lower triangle of the distance matrix for all unique
	# combinations of two strains
	for (i in 2:nrow(res)) {
		for (j in 1:(i-1)) {
			seqs_clean <- remove_ambiguous_residues(
				c(seqs[[i]], seqs[[j]]),
				ambiguous_residues
			)
			res[[i, j]] <- pepitope(
				seqs_clean[[1]], seqs_clean[[2]],
				subtype, mode, harmonize_b_lineages, ambiguous_residues
			)
		}
	}
	
	# Set the diagonal to zero
	diag(res) <- 0
	
	# Quickly fill in the upper triangle
	out <- res |> as.dist() |> as.matrix()
	
	return(out)
}

# From substitution metrics script
# https://github.com/ahgroup/influenza-antigenic-distance-data/blob/33d9b7607176207cc37de23df1c7bfa84705673b/R/substitution.R

# This function generates an R matrix containing the pairwise Grantham
# distances from the table on wikipedia.
# TODO cite the original paper(s)
# Note that this function gets run every time the Grantham distance is
# calculated. Because this function is very fast, that's probably faster than
# saving this to a file and loading it every time, which depends on the I/O
# speed of the drive.
generate_grantham_matrix <- function() {
	wiki_aa_order <- c(
		"s", "r", "l", "p", "t", "a", "v", "g", "i", "f", "y",
		"c", "h", "q", "n", "k", "d", "e", "m", "w"
	)
	
	grantham_matrix <-
		matrix(
			ncol = 20,
			nrow = 20,
			dimnames = list(wiki_aa_order, wiki_aa_order)
		)
	grantham_values <-
		c(
			110, 145, 74, 58, 99, 124, 56, 142, 155, 144, 112, 89, 68, 46, 121, 65,
			80, 135, 177, 102, 103, 71, 112, 96, 125, 97, 97, 77, 180, 29, 43,
			86, 26, 96, 54, 91, 101, 98, 92, 96, 32, 138, 5, 22, 36, 198, 99, 113,
			153, 107, 172, 138, 15, 61, 38, 27, 68, 42, 95, 114, 110, 169, 77, 76, 
			91, 103, 108, 93, 87, 147, 58, 69, 59, 89, 103, 92, 149, 47, 42, 65, 78,
			85, 65, 81, 128, 64, 60, 94, 113, 112, 195, 86, 91, 111, 106, 126, 107,
			84, 148, 109, 29, 50, 55, 192, 84, 96, 133, 97, 152, 121, 21, 88, 135,
			153, 147, 159, 98, 87, 80, 127, 94, 98, 127, 184, 21, 33, 198, 94, 109,
			149, 102, 168, 134, 10, 61, 22, 205, 100, 116, 158, 102, 177, 140, 28,
			40, 194, 83, 99, 143, 85, 160, 122, 36, 37, 174, 154, 139, 202, 154, 170,
			196, 215, 24, 68, 32, 81, 40, 87, 115, 46, 53, 61, 29, 101, 130, 94,
			23, 42, 142, 174, 101, 56, 95, 110, 45, 160, 181, 126, 152, 67
		)
	# R will fill in the upper triangle by COLUMNS instead of by ROWS, which
	# would be wrong with the order the values are typed in the vector above.
	# Fix this by filling in the lower triangle instead and then transposing
	# which gives the correct order.
	grantham_matrix[lower.tri(grantham_matrix, diag = FALSE)] <- grantham_values
	grantham_matrix <- t(grantham_matrix)
	diag(grantham_matrix) <- rep(0, times = 20)
	grantham_matrix[lower.tri(grantham_matrix)] <-
		t(grantham_matrix)[lower.tri(grantham_matrix)]
	
	return(grantham_matrix)
}

# Function to generate the FLU substitution matrix
generate_FLU_matrix <- function() {
	# The input data are downloaded from the linked provided in
	# Dang, C.C., Le, Q.S., Gascuel, O. et al. FLU, an amino acid substitution
	# model for influenza proteins. BMC Evol Biol 10, 99 (2010).
	# https://doi.org/10.1186/1471-2148-10-99
	# Which is
	# ftp://ftp.sanger.ac.uk/pub/1000genomes/lsq/FLU
	# This is copied and pasted from the provided file and formatted for R
	# with no other changes other than adding commas
	# So we don't have to read/write every time
	# THe provided file is included in the data/raw directory as
	# "flu-matrix-PAML.txt" for anyone who wants to verify that these match.
	input_data <- c(
		0.138658764751059,
		0.0533665787145181, 0.161000889039552,
		0.584852305649886, 0.00677184253227681, 7.73739287051356,
		0.0264470951166826, 0.16720700818221, 1.30249856764315e-005, 0.014132062548787,
		0.353753981649393, 3.29271694159791, 0.530642655337477, 0.145469388422239, 0.00254733397966779,
		1.4842345032161, 0.124897616909194, 0.0616521921873234, 5.37051127867923, 3.91106992668137e-011, 1.19562912226203,
		1.13231312248046, 1.19062446519178, 0.322524647863997, 1.93483278448943, 0.116941459124876, 0.108051341246072, 1.59309882471598,
		0.214757862168721, 1.87956993845887, 1.38709603234116, 0.887570549414031, 0.0218446166959521, 5.33031341222104, 0.256491863423002, 0.0587745274250666,
		0.149926734229061, 0.246117171830255, 0.21857197541607, 0.0140859174993809, 0.00111215807314139, 0.0288399502994541, 0.0142107118685268, 1.62662283098296e-005, 0.243190142026506,
		0.0231169515264061, 0.296045557460629, 0.000835873174542931, 0.00573068208525287, 0.00561362724916376, 1.02036695531654, 0.016499535540562, 0.00651622937676521, 0.321611693603646, 3.51207228207807,
		0.474333610192982, 15.3000966197798, 2.6468479652886, 0.290042980143818, 3.83228119049152e-006, 2.559587177122, 3.88148880863814, 0.264148929349066, 0.347302791211758, 0.227707997165566, 0.129223639195248,
		0.0587454231508643, 0.890162345593224, 0.00525168778853117, 0.0417629637305017, 0.111457310321926, 0.190259181297527, 0.313974351356074, 0.00150046692269255, 0.00127350890508147, 9.01795420287895, 6.74693648486614, 1.33129161941264,
		0.0804909094320368, 0.0160550314767596, 0.000836445615590923, 1.0600102849456e-006, 0.10405366623526, 0.0326806570137471, 0.00100350082518749, 0.00123664495412902, 0.119028506158521, 1.46335727834648, 2.98680003596399, 0.319895904499071, 0.279910508981581,
		0.659311477863896, 0.154027179890711, 0.0364417719063219, 0.188539456415654, 1.59312060172652e-013, 0.712769599068934, 0.319558828428154, 0.0386317614553493, 0.924466914225534, 0.0805433268150369, 0.634308520867322, 0.195750631825315, 0.0568693216513547, 0.0071324304661639,
		3.01134451903854, 0.950138410087378, 3.88131053061457, 0.338372183381345, 0.336263344504404, 0.487822498528951, 0.307140298031341, 1.58564657669139, 0.580704249811294, 0.290381075260226, 0.570766693213698, 0.283807671568883, 0.00702658828739369, 0.996685669575839, 2.08738534433198,
		5.4182981753166, 0.183076905018197, 2.14033231636063, 0.135481232622983, 0.011975265782196, 0.60234096342392, 0.2801248951174, 0.0188080299490973, 0.368713573381758, 2.90405228596936, 0.0449263566753846, 1.52696419998775, 2.03151132062208, 0.000134906239484254, 0.54225109402693, 2.2068599339404,
		0.195966354027106, 1.36942940801512, 0.000536284040016542, 1.4893873721753e-005, 0.0941066800969967, 0.0440205200833047, 0.155245492137294, 0.196486447133033, 0.0223729191088972, 0.0321321499585514, 0.431277662888057, 4.97641445484395e-005, 0.0704600385245663, 0.814753093809928, 0.000431020702277328, 0.0998357527014247, 0.207066205546908,
		0.0182892882245349, 0.0998554972524385, 0.373101926513925, 0.525398542949365, 0.601692431136271, 0.0722059354079545, 0.104092870343653, 0.0748149970972622, 6.44895444648517, 0.273934263183281, 0.340058468374384, 0.0124162215506117, 0.874272174533394, 5.39392424532822, 0.000182294881489116, 0.392552239890831, 0.124898020409882, 0.42775543040588,
		3.53200526987468, 0.103964386383736, 0.0102575172450253, 0.297123975243582, 0.0549045639492389, 0.406697814049488, 0.285047948309311, 0.337229618868315, 0.0986313546653266, 14.3940521944257, 0.890598579382591, 0.0731279296372675, 4.90484223478739, 0.592587985458668, 0.0589719751511691, 0.0882564232979724, 0.654109108255219, 0.256900461407996, 0.167581646770807,
		0.0470718, 0.0509102, 0.0742143, 0.0478596, 0.0250216, 0.0333036, 0.0545874, 0.0763734, 0.0199642, 0.0671336, 0.0714981, 0.0567845, 0.0181507, 0.0304961, 0.0506561, 0.0884091, 0.0743386, 0.0185237, 0.0314741, 0.0632292
	)
	
	n_entries <- length(input_data)
	
	# The data are in PAML/phyML format, which means that the data are formatted
	# as first the lower triangle (no diagonal) of the substitution matrix,
	# followed by a blank line, followed by the list of equilibrium frequencies.
	# So first we can remove the last 20 elements of the input, which are the
	# equilibrium frequencies.
	equilibrium_frequencies <- tail(input_data, 20)
	matrix_data <- head(input_data, n_entries - 20)
	
	# PAML format lists the amino acids in alphabetical order
	paml_aa_order <- c(
		"a", "r", "n", "d", "c", "q", "e", "g", "h", "i", "l",
		"k", "m", "f", "p", "s", "t", "w", "y", "v"
	)
	
	# Create an empty matrix and fill in the upper triangle with the matrix
	# entries, then transpose the matrix. We have to do it this way instead of
	# just filling in the lower triangle cause this ensures the matrix is filled
	# in the correct order (by row instead of by column)
	sub_mat <- matrix(
		nrow = 20, ncol = 20,
		dimnames = list(paml_aa_order, paml_aa_order)
	)
	sub_mat[upper.tri(sub_mat)] <- matrix_data
	out <- t(sub_mat)
	
	# Now fill in the diagonal and upper triangle quickly
	out <- as.matrix(as.dist(out))
	
	return(out)
}

# Calculate the distance between two amino acid strings based on a substitution
# matrix. Right now this only supports Grantham's distance.
substitution <- function(seq1, seq2, method = "grantham", ambiguous_residues = "xX?-") {
	# Validate input sequences
	if (!is.character(seq1) || !is.character(seq2)) {
		stop("Both seq1 and seq2 must be character vectors.")
	}
	
	# Make sure sequences are the same length
	if (length(seq1) != length(seq2)) {
		stop("Your sequences should be aligned! They are different lengths.")
	}
	
	# Remove ambiguous residues from the two sequences
	seqs <- c(seq1, seq2) |>
		tolower() |>
		remove_ambiguous_residues(ambiguous_residues)
	seq1_ua <- strsplit(seqs[[1]], "")[[1]]
	seq2_ua <- strsplit(seqs[[2]], "")[[1]]
	
	# Sanitize method argument
	method_san <- tolower(method)
	
	# Load the correct substitution matrix
	substitution_matrix <- switch(
		method_san,
		"grantham" = generate_grantham_matrix(),
		"flu" = generate_FLU_matrix(),
		stop(paste0(
			"Invalid 'method' specified. You specified: '", method, "', but the ",
			"allowed options are:\n'Grantham', 'FLU' (case-insensitive)."
		))
	)
	
	# Ensure sequences contain valid amino acids
	valid_aa <- rownames(substitution_matrix)
	all_seq_chars <- union(seq1_ua, seq2_ua)
	unallowed_chars <- setdiff(all_seq_chars, valid_aa)
	
	if (length(unallowed_chars) > 0) {
		stop(paste0(
			"Sequences must only contain valid amino acid residues.",
			" One of your sequences has the prohibited character(s): ",
			unallowed_chars
		))
	}
	
	# Compute the distances from the substitution matrix, use the mean so
	# it's normalized by protein length.
	substitution_distance <- substitution_matrix[cbind(seq1_ua, seq2_ua)]
	
	# Return the total substitution distance normalized by the sequence lengths.
	total_distance <- mean(substitution_distance)
	
	return(total_distance)
}

# Now apply to all pairs in a vector of sequences
# TODO this seems like it could be generalized to accept a distance function
# as an argument, see e.g. p-epitope.
# We do it this way with the loop because even though outer() should probably
# be faster this only does n * (n-1) / 2 calculations instead of n ^ 2.
# Not sure which is faster.
dist_substitution <- function(seqs,  method = "grantham", ambiguous_residues = "xX?") {
	# Set up an empty matrix to hold results
	res <- matrix(
		nrow = length(seqs),
		ncol = length(seqs),
		dimnames = list(names(seqs), names(seqs))
	)
	
	# Calculate the lower triangle of the distance matrix for all unique
	# combinations of two strains
	for (i in 2:nrow(res)) {
		for (j in 1:(i-1)) {
			res[[i, j]] <- substitution(
				seqs[[i]], seqs[[j]],
				method, ambiguous_residues
			)
		}
	}
	
	# Set the diagonal to zero
	diag(res) <- 0
	
	# Quickly fill in the upper triangle
	out <- res |> as.dist() |> as.matrix()
	
	return(out)
	
	return(dist)
}

# Cleaning function
# https://github.com/ahgroup/influenza-antigenic-distance-data/blob/33d9b7607176207cc37de23df1c7bfa84705673b/R/distance-calculation.R

tidy_joined_distances <- function(combined_distance_data) {
	tidy_distance_data <-
		combined_distance_data |>
		# First pivot longer so that instead of one list-column for each distance
		# matrix, we have one column specifying the metric and another list-column
		# containing all of the distance matrices.
		tidyr::pivot_longer(
			cols = dplyr::contains("dist_"),
			names_to = "metric",
			values_to = "distance_df"
		) |>
		# Now we can do some cleanup of the ID variables (subtype and metric),
		# and we map the tidy_dist_mat() function over the list-column of distance
		# matrices. This converts each of the distance matrices into tidy format
		# which only has three columns.
		dplyr::mutate(
			type_subtype = factor(type_subtype),
			metric = factor(gsub("dist_", "", metric)),
			distance_df = purrr::map(distance_df, tidy_dist_mat)
		) |>
		# Now we can finally unnest the distance matrix column, which means that
		# instead of having a list-column of distance data frames, we will have one
		# rectangular data frame which has all of the distances in it.
		tidyr::unnest(distance_df) |>
		# Make the missing cartographic data explicit -- this makes the d values
		# new rows with NA values for combinations that don't exist in the data.
		# The nesting() part ensures that only subtype/Strain1/Strain2 combinations
		# that already exist in the data are used, so we don't get nonsensical
		# stuff like mismatched strains and subtypes.
		tidyr::complete(tidyr::nesting(type_subtype, Strain1, Strain2), metric)
	
	return(tidy_distance_data)
}

minmax <- function(x) {
	return((x - min(x)) / (max(x) - min(x)))
}
