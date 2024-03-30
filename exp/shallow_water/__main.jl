import Dates


function Shallow_Water_Main(;model_name::String = "Shallow_Water",
                            radius::Float64 = 6371.0e3,
                            omega::Float64 = 7.292e-5,
                            grav::Float64 = 9.80,
                            sea_level_ps_ref::Float64 = 1.0e5,
                            hour_to_sec::Int64 = 3600,
                            day_to_sec::Int64 = 86400,
                            init_step::Bool = true)
    """
    TODO::
    """
    ########################################
    ### Model setting
    
    # Meta
    creation_time = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    fname_prefix = "SWM_"
    fname = string(rand(1000000:10000000-1))
    fname_suffix = ".jld2" # JLD2 package
    output_filename = fname_prefix * fname * fname_suffix 
    
    # Resolution
    nλ = 512
    nθ = 256
    nd = 1
    num_fourier = floor(Int64, nθ*(2/3))
    num_spherical = num_fourier + 1
    
    # Hyper-viscosity
    damping_order = 4
    damping_coef = 1.e-05
    robert_coef  = 0.04 
    implicit_coef = 0.5
    
    # Time
    start_time = 0
    end_time = 86400*5
    Δt = 60
    
    # Convective efficiency
    kappa_chi = 1 / (0.5 * 86400)
    kappa_lambda = 900.0 * grav
    
    # Gravity wave speed
    hbar_grav_wave = 3000.0 * grav
    
    # Background
    vor_b_amplitude = 6.2e-5
    vor_b_latitude = 19.0
    vor_b_width = 3.0
    
    # Perturbation
    vor_a_amplitude = 6e-5
    vor_a_longitude = 180.0
    vor_a_latitude = 20.0
    vor_a_width = 2.5
    
    ########################################
    ### Initialize model
    
    # spectral_spherical_mesh
    mesh = Spectral_Spherical_Mesh(num_fourier = num_fourier, 
                                   num_spherical = num_spherical,
                                   nλ = nλ, 
                                   nθ = nθ, 
                                   nd = nd,
                                   radius = radius)
    
    # vertical_coordinate
    θc, λc = mesh.θc,  mesh.λc
    cosθ, sinθ = mesh.cosθ, mesh.sinθ
    vert_coord = Vert_Coordinate(nλ = nλ, 
                                 nθ = nθ, 
                                 nd = nd,
                                 vert_coord_option = "even_sigma", 
                                 vert_difference_option = "simmons_and_burridge", 
                                 vert_advect_scheme = "second_centered_wts",
                                 p_ref = sea_level_ps_ref)
    
    # atmospheric data
    atmo_data = Atmo_Data(name = model_name, 
                          nλ = nλ, 
                          nθ = nθ, 
                          nd = nd,
                          do_mass_correction = false,
                          do_energy_correction = false,
                          do_water_correction = false,
                          use_virtual_temperature = false,
                          sinθ = sinθ,
                          radius = radius, 
                          omega = omega,
                          grav = grav)
    
    # dynamical data
    dyn_data = Dyn_Data(name = model_name,
                        num_fourier = num_fourier, 
                        num_spherical = num_spherical,
                        nλ = nλ, 
                        nθ = nθ, 
                        nd = nd)
    
    # integrator
    integrator = Filtered_Leapfrog(robert_coef = robert_coef, 
                                   damping_order = damping_order, 
                                   damping_coef = damping_coef,
                                   eigen = mesh.laplacian_eig, 
                                   implicit_coef = implicit_coef,
                                   Δt = Δt, 
                                   init_step = init_step, 
                                   start_time = start_time, 
                                   end_time = end_time)
    
    # output manager
    output_manager = Output_Manager(mesh = mesh, 
                                    integrator = integrator, 
                                    vert_coord = vert_coord, 
                                    start_time = start_time, 
                                    end_time = end_time, 
                                    creation_time = creation_time, 
                                    filename = output_filename)
    
    ########################################
    ### Initialize variable field

    # Main field
    grid_u = dyn_data.grid_u_c 
    grid_v = dyn_data.grid_v_c 
    grid_vor = dyn_data.grid_vor 
    grid_div = dyn_data.grid_div 
    grid_h = dyn_data.grid_ps_c
    # grid_q = dyn_data.grid_q_c
    spe_vor_c = dyn_data.spe_vor_c
    spe_div_c = dyn_data.spe_div_c
    spe_h_c = dyn_data.spe_lnps_c
    # spe_q_c = dyn_data.spe_q_c
    
    # Background
    background_field = Background_Vorticity_Strip(mesh = mesh,
                                                  atmo_data = atmo_data,
                                                  dyn_data = dyn_data,
                                                  vor_amp = vor_b_amplitude,
                                                  vor_lat = vor_b_latitude, 
                                                  vor_width = vor_b_width,
                                                  IsDisplay = true)
    # Perturbation
    perturbation_field = Isolated_Vorticity_Blob(mesh = mesh,
                                                 atmo_data = atmo_data,
                                                 dyn_data = dyn_data,
                                                 vor_amp = vor_a_amplitude,
                                                 vor_lon = vor_a_longitude,
                                                 vor_lat = vor_a_latitude,
                                                 vor_width = vor_a_width,
                                                 IsDisplay = true)
    
    # Total
    grid_u .= (background_field[1] + perturbation_field[1])
    grid_v .= (background_field[2] + perturbation_field[2]) 
    grid_vor .= (background_field[3] + perturbation_field[3])
    grid_div .= (background_field[4] + perturbation_field[4])
    grid_h .= (background_field[5] + perturbation_field[5]) .+ hbar_grav_wave
    # grid_q = 0.0
    spe_vor_c .= (background_field[6] + perturbation_field[6])
    spe_div_c .= (background_field[7] + perturbation_field[7])
    spe_h_c .= (background_field[8] + perturbation_field[8])
    spe_h_c[1,1] += hbar_grav_wave 
    # spe_q_c .= 0.0
    
    #
    Update_Output_Init!(output_manager = output_manager, 
                        dyn_data = dyn_data, 
                        current_time = integrator.time)
    
    ########################################
    ### Simulation stage
    while (integrator.time+Δt) <= end_time
        # 1
        Shallow_Water_Physics!(dyn_data = dyn_data, 
                               kappa_chi = kappa_chi, 
                               kappa_lambda = kappa_lambda)
        Shallow_Water_Dynamics!(mesh = mesh, 
                                atmo_data = atmo_data, 
                                dyn_data = dyn_data, 
                                integrator = integrator, 
                                h_0 = hbar_grav_wave)
        # 2
        Time_Advance!(dyn_data = dyn_data)
        integrator.time += Δt
        # 3
        Update_Output!(output_manager = output_manager, 
                       dyn_data = dyn_data, 
                       current_time = integrator.time)
        #
        if integrator.init_step
            Update_Init_Step!(integrator = integrator)
        end
        
        if (integrator.time%hour_to_sec == 0)
            println(repeat("***", 30))
            println("Hour: ", (÷(integrator.time, hour_to_sec)))
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
        end
    end
    
    ########################################
    ### Output stage
    Finalize_Output!(output_manager = output_manager)
    Generate_Output!(output_manager = output_manager)
    
    return nothing
end