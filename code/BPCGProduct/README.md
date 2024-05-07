To run the experiments:

1. Set up experimental parameters in `config.yml`
2. decide LMOs in `run_experiments.jl`
3. Run 
    ```bash
    julia --project=. test/run_experiments.jl > <your_output_logfile>.log 2>&1
    ```
    from terminal to have the experiments' log printed to custom logfile