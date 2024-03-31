export Display_Initial_Background!
export Display_Initial_Perturbation!
export Display_Current_State!


function Display_Initial_Background!(;logpath::String,
                                     grid_u::Array{Float64,3},
                                     grid_v::Array{Float64,3},
                                     grid_h::Array{Float64,3},
                                     grid_vor::Array{Float64,3},
                                     grid_div::Array{Float64,3})
    
    # Display on terminal
    println(repeat("###", 30))
    println("Initial background      zonal wind: ",
            (round(minimum(grid_u); digits = 4), 
            round(maximum(grid_u); digits = 4)),
            " (m/s)")
    println("Initial background meridional wind: ",
            (round(minimum(grid_v); digits = 4), 
            round(maximum(grid_v); digits = 4)),
            " (m/s)")
    println("Initial background   geopot height: ", 
            (round(minimum(grid_h); digits = 4), 
            round(maximum(grid_h); digits = 4)),
            " (m * m*s-2)")
    println("Initial background       vorticity: ", 
            (round(minimum(grid_vor); digits = 9), 
            round(maximum(grid_vor); digits = 9)),
            " (1/s)")
    println("Initial background      divergence: ", 
            (round(minimum(grid_div); digits = 9), 
            round(maximum(grid_div); digits = 9)),
            " (1/s)")
    println(repeat("###", 30))
    
    # Display on file
    open(logpath, "a+") do file
        #
        write(file, repeat("###", 30), "\n")
        #
        write(file, "Initial background      zonal wind: ")
        write(file, "( ", string(round(minimum(grid_u); digits = 4)), ", ")
        write(file, string(round(maximum(grid_u); digits = 4)), " )")
        write(file, " (m/s)\n")
        #
        write(file, "Initial background meridional wind: ")
        write(file, "( ", string(round(minimum(grid_v); digits = 4)), ", ")
        write(file, string(round(maximum(grid_v); digits = 4)), " )")
        write(file, " (m/s)\n")
        #
        write(file, "Initial background   geopot height: ")
        write(file, "( ", string(round(minimum(grid_h); digits = 4)), ", ")
        write(file, string(round(maximum(grid_h); digits = 4)), " )")
        write(file, " (m * m*s-2)\n")
        #
        write(file, "Initial background       vorticity: ")
        write(file, "( ", string(round(minimum(grid_vor); digits = 9)), ", ")
        write(file, string(round(maximum(grid_vor); digits = 9)), " )")
        write(file, " (1/s)\n")
        #
        write(file, "Initial background      divergence: ")
        write(file, "( ", string(round(minimum(grid_div); digits = 9)), ", ")
        write(file, string(round(maximum(grid_div); digits = 9)), " )")
        write(file, " (1/s)\n")
        #
        write(file, repeat("###", 30), "\n")
    end
end


function Display_Initial_Perturbation!(;logpath::String,
                                       grid_u::Array{Float64,3},
                                       grid_v::Array{Float64,3},
                                       grid_h::Array{Float64,3},
                                       grid_vor::Array{Float64,3},
                                       grid_div::Array{Float64,3})
    
    # Display on terminal
    println(repeat("###", 30))
    println("Initial perturbation      zonal wind: ",
            (round(minimum(grid_u); digits = 4), 
            round(maximum(grid_u); digits = 4)),
            " (m/s)")
    println("Initial perturbation meridional wind: ",
            (round(minimum(grid_v); digits = 4), 
            round(maximum(grid_v); digits = 4)),
            " (m/s)")
    println("Initial perturbation   geopot height: ", 
            (round(minimum(grid_h); digits = 4), 
            round(maximum(grid_h); digits = 4)),
            " (m * m*s-2)")
    println("Initial perturbation       vorticity: ", 
            (round(minimum(grid_vor); digits = 9), 
            round(maximum(grid_vor); digits = 9)),
            " (1/s)")
    println("Initial perturbation      divergence: ", 
            (round(minimum(grid_div); digits = 9), 
            round(maximum(grid_div); digits = 9)),
            " (1/s)")
    println(repeat("###", 30))
    
    # Display on file
    open(logpath, "a+") do file
        #
        write(file, repeat("###", 30), "\n")
        #
        write(file, "Initial perturbation      zonal wind: ")
        write(file, "( ", string(round(minimum(grid_u); digits = 4)), ", ")
        write(file, string(round(maximum(grid_u); digits = 4)), " )")
        write(file, " (m/s)\n")
        #
        write(file, "Initial perturbation meridional wind: ")
        write(file, "( ", string(round(minimum(grid_v); digits = 4)), ", ")
        write(file, string(round(maximum(grid_v); digits = 4)), " )")
        write(file, " (m/s)\n")
        #
        write(file, "Initial perturbation   geopot height: ")
        write(file, "( ", string(round(minimum(grid_h); digits = 4)), ", ")
        write(file, string(round(maximum(grid_h); digits = 4)), " )")
        write(file, " (m * m*s-2)\n")
        #
        write(file, "Initial perturbation       vorticity: ")
        write(file, "( ", string(round(minimum(grid_vor); digits = 9)), ", ")
        write(file, string(round(maximum(grid_vor); digits = 9)), " )")
        write(file, " (1/s)\n")
        #
        write(file, "Initial perturbation      divergence: ")
        write(file, "( ", string(round(minimum(grid_div); digits = 9)), ", ")
        write(file, string(round(maximum(grid_div); digits = 9)), " )")
        write(file, " (1/s)\n")
        #
        write(file, repeat("###", 30), "\n") 
    end
end


function Display_Current_State!(;logpath::String,
                                tick::Int64,
                                grid_u::Array{Float64,3},
                                grid_v::Array{Float64,3},
                                grid_h::Array{Float64,3},
                                grid_vor::Array{Float64,3},
                                grid_div::Array{Float64,3})
    # Display on terminal
    println(repeat("***", 30))
    println("Hour: ", string(tick))
    println("     zonal wind: ",
            (round(minimum(grid_u); digits = 4), 
            round(maximum(grid_u); digits = 4)),
            " (m/s)")
    println("meridional wind: ",
            (round(minimum(grid_v); digits = 4), 
            round(maximum(grid_v); digits = 4)),
            " (m/s)")
    println("  geopot height: ", 
            (round(minimum(grid_h); digits = 4), 
            round(maximum(grid_h); digits = 4)),
            " (m * m*s-2)")
    println("      vorticity: ", 
            (round(minimum(grid_vor); digits = 9), 
            round(maximum(grid_vor); digits = 9)),
            " (1/s)")
    println("     divergence: ", 
            (round(minimum(grid_div); digits = 9), 
            round(maximum(grid_div); digits = 9)),
            " (1/s)")
    println(repeat("***", 30))
    
    # Display on file
    open(logpath, "a+") do file
        #
        write(file, repeat("***", 30), "\n")
        #
        write(file, "Hour: ", string(tick), "\n")
        #
        write(file, "     zonal wind: ")
        write(file, "( ", string(round(minimum(grid_u); digits = 4)), ", ")
        write(file, string(round(maximum(grid_u); digits = 4)), " )")
        write(file, " (m/s)\n")
        #
        write(file, "meridional wind: ")
        write(file, "( ", string(round(minimum(grid_v); digits = 4)), ", ")
        write(file, string(round(maximum(grid_v); digits = 4)), " )")
        write(file, " (m/s)\n")
        #
        write(file, "  geopot height: ")
        write(file, "( ", string(round(minimum(grid_h); digits = 4)), ", ")
        write(file, string(round(maximum(grid_h); digits = 4)), " )")
        write(file, " (m * m*s-2)\n")
        #
        write(file, "      vorticity: ")
        write(file, "( ", string(round(minimum(grid_vor); digits = 9)), ", ")
        write(file, string(round(maximum(grid_vor); digits = 9)), " )")
        write(file, " (1/s)\n")
        #
        write(file, "     divergence: ")
        write(file, "( ", string(round(minimum(grid_div); digits = 9)), ", ")
        write(file, string(round(maximum(grid_div); digits = 9)), " )")
        write(file, " (1/s)\n")
        write(file, repeat("***", 30), "\n")
    end
end