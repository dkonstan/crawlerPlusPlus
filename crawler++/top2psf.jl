using ArgParse


function top2psf(top, psf)
	# converts crawler.jl topology to a minimal psf topology file that VMD can read

	bondIdx = []
	open(top, "r") do top
		contents = readlines(top)
		i = 1
		while i <= length(contents)
			if occursin("<bonds>", contents[i])
				i += 1
				while !occursin("<end>", contents[i])
					bondContents = split(contents[i], " ")
					bond1 = parse(Int64, bondContents[1])
					bond2 = parse(Int64, bondContents[2])
					bondLength = parse(Float64, bondContents[3])
					bondK = parse(Float64, bondContents[4])
					push!(bondIdx, [bond1, bond2])
					i += 1
				end
			end
			i += 1
		end

	end

	nAtoms = maximum(hcat(bondIdx...)') + 1  # + 1 for crawler++ because starts with 0
	nBonds = size(bondIdx, 1)
	open(psf, "w+") do top
		write(top, "PSF\n")
		write(top, "     0 !NTITLE\n")
		write(top, "     $(nAtoms) !NATOM\n")
		for i in 1:nAtoms
			write(top, "       $(i) MAIN 1    X    C    C      0.000000       0.00000           0\n")
		end
		write(top, "\n")
		write(top, "     $(nBonds) !NBOND: bonds\n")
		write(top, "       ")
		for i in 1:nBonds
			# 1 and 2 digit numbers only for now
			write(top, "$(bondIdx[i][1] + 1)")
			write(top, repeat(" ", 7 - length("$(bondIdx[i][1])")))
			write(top, "$(bondIdx[i][2] + 1)")
			write(top, repeat(" ", 7 - length("$(bondIdx[i][2])")))
			if i % 4 == 0
				write(top, "\n       ")
			end
		end
	end
end


argTable = ArgParseSettings()

@add_arg_table argTable begin
    "topology"
        arg_type = String
        required = true
    "psf"
        arg_type = String
        required = true
  	end

args = parse_args(argTable)

top2psf(args["topology"], args["psf"])
