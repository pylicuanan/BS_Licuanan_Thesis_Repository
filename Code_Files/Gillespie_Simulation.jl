using DifferentialEquations, Plots, LinearAlgebra, Statistics, Catalyst, CSV, DataFrames

default(; lw = 2)

# Define the reaction network
NCAX_model = @reaction_network begin
    p_bar, N + C --> 2C
    r, C --> A
    u, A --> X
    m, N + X --> 2N
    c_bar, C + X --> 2C
end 

Final_Time = 200.0
Total_Population = 400
Total_Iteration = 400
saveat_interval = 1.0

# Define parameter ranges
P = (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0)
Rbar = (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0)

u_values = 1


Repetition = range(301, Total_Iteration, step = 1)

Final_N_mean = []
Final_C_mean = []
Final_A_mean = []
Final_X_mean = []

for x in P, y in Rbar, z in u_values
    N_values = [] 
    C_values = []
    A_values = []
    X_values = []
    
    for i in Repetition        
        p = (:p_bar => x / Total_Population, :r => y, :u => z, :m => 0.5 / Total_Population, :c_bar => ((0.5)(2) / Total_Population))
        u₀ = [:N => 399, :C => 1, :A => 0, :X => 0]
        tspan = (0.0, Final_Time)
        prob = DiscreteProblem(NCAX_model, u₀, tspan, p)

        jump_prob = JumpProblem(NCAX_model, prob, Direct(); save_positions = (false, false))
        
        # Solve the problem
        sol = solve(jump_prob, SSAStepper(); saveat = saveat_interval)
        ts = sol.t
        N_vals = getindex.(sol.u, 1) ./ Total_Population
        C_vals = getindex.(sol.u, 2) ./ Total_Population
        A_vals = getindex.(sol.u, 3) ./ Total_Population
        X_vals = getindex.(sol.u, 4) ./ Total_Population
        
        # Plot the solution
        display(plot(ts, [N_vals C_vals A_vals X_vals], title="Population Plot rep = $i, p = $x, r = $y, u = $z", 
                        xlabel="Timestep", ylabel="Population", 
                        label=["N(t)" "C(t)" "A(t)" "X(t)"],
                        color=[:green :red :blue :gray]))
        
        filename = "C:/Users/phili/OneDrive/Documents/Research/Julia/V4_NCAX_Repetition/popplot_rep=$i, p=$x, r=$y, u=$z.png"
        savefig(filename)

        # Collect individual values
        N_value = map(x -> x[1], sol.u)
        push!(N_values, N_value)
        
        C_value = map(x -> x[2], sol.u)
        push!(C_values, C_value)
        
        A_value = map(x -> x[3], sol.u)
        push!(A_values, A_value)
        
        X_value = map(x -> x[4], sol.u) 
        push!(X_values, X_value)



        # Save individual population data
        individual_data = DataFrame(Time=collect(0:saveat_interval:Final_Time), N=N_value, C=C_value, A=A_value, X=X_value)
        csv_filename = "C:/Users/phili/OneDrive/Documents/Research/Julia/V4_CSVData/individual_populations_p=$x, r=$y, u=$z, rep=$i.csv"
        CSV.write(csv_filename, individual_data)
    end
end


