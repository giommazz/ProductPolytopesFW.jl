# `slurm_loop.jl`
#!/usr/bin/env julia
using ProductPolytopesAFW
using Printf

# ***********************************************************************
# HELPERS FUNCTIONS
# ***********************************************************************

function write_modified_config(config::Config, k::Integer, n::Integer, seed::Integer, dir::AbstractString)
    modified_config = modify_config(config, k=k, n=n, seed=seed)
    config_filename = joinpath(dir, "config_k$(k)_n$(n)_s$(seed).yml")
    return write_config(modified_config, config_filename)
end

# replaces the line "config = Config("examples/config.yml")" in `script_filename` with a new line "config = Config(config_filename)" and saves text into new file
function new_modified_file(
    filename::AbstractString,
    src_textline::AbstractString,
    replacement_textline::AbstractString,
    k::Integer,
    n::Integer,
    seed::Integer)

    # get direcgtory name from `filename`
    dir = dirname(filename)
    
    # read the template `script_filename`
    src_text = read(filename, String)

    # replace the `src_textline` with `replacement_textline`
    clean  = x -> strip(replace(x, "\r\n" => "\n"))
    new_text = replace(src_text, clean(src_textline) => replacement_textline)

    # craft output path and write new file
    # returns name of `script_filename` and format
    scriptname, format  = splitext(basename(filename))[1], splitext(basename(filename))[2]
    new_filename    = joinpath(dir, "$(scriptname)_k$(k)_n$(n)_s$(seed)$(format)")
    mkpath(dirname(new_filename))  # be sure the dir exists
    write(new_filename, new_text)

    return abspath(new_filename)
end

"""
    safe_submit(cmd::Cmd)

Runs `cmd` asynchronously (`wait=false`) on SLURM.  If `sbatch` errors, `safe_submit` catches the exception and only emits a warning
"""
function safe_submit(cmd::Cmd)
    try
        # `wait=false`: Julia will launch the external command `cmd` and then give control back to `safe_submit`
        #       without sitting around waiting for that command to finish.
        #       By default, `run(cmd)` blocks until the external program exits, but adding `wait=false` makes it asynchronous
        #       Enclosing everything in a try-catch statement ensures any SLURM problem will be caught without blockng the Julia script
        #       So, with this, SLURM enqueues the job and Julia proceeds to the next step
        run(cmd)
    catch err
        @warn "slurm submission failed" cmd=cmd error=err
    end
end

function main(dir, config, seed)
    for k in list_k, n in list_n
        
        # write new modified config "examples/config.yml", with modified parameters k, n, seed
        new_config_filename = write_modified_config(config, k, n, seed, dir)

        # write new modified script "examples/compute_intersection_custom_full_warmup_slurm.jl", with new line "config = Config($(new_config_filename))"
        new_jl_experimentsscript_filename = new_modified_file(
            joinpath(dir, "compute_intersection_custom_full_warmup.jl"), 
            "config = Config(\"examples/config.yml\")",
            "config = Config(\"$(new_config_filename)\")",
            k, n, seed)
        println(new_jl_experimentsscript_filename)

        # write new modified script "examples/slurm_experiments.sh", with new line "config = Config($(new_config_filename))"
        new_sh_slurmscript_filename = new_modified_file(
            joinpath(dir, "slurm_experiments.sh"),
            "#SBATCH --cpus-per-task=2       # reserve a certain amount of cores for this job",
            "#SBATCH --cpus-per-task=$(k)       # reserve a certain amount of cores for this job",
            k, n, seed)
        println(new_sh_slurmscript_filename)
        
        # run `sbatch` command
        slurm_command = `sbatch $(new_sh_slurmscript_filename) $(new_jl_experimentsscript_filename) $dir/results_linesearch_afw/ $new_config_filename`
        println("Submitting job with:\n  ", slurm_command)
        # launch command, catch any errors
        # safe_submit(slurm_command)
        println()

        seed += 1
    end
end



# ***********************************************************************
# INSTANCE PARAMETERS FOR THE RUNS
# ***********************************************************************
list_k    = [5]           # number of polytopes
list_n    = [10000, 10000]      # dimension of each polytope
seed = 344               # starting seed, will be incremented in the loop
dir = "examples"
config = Config(joinpath(dir, "config.yml"))

main(dir, config, seed)
