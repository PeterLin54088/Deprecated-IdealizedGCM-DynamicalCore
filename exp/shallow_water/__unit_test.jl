using LinearAlgebra
# derivative will lose accuracy, which depends on the last two spherical modes


function IC_Profiler(IC::String = "Case 1")
    """
    TODO::Thickness IC might needed
    """
    # Earth
    radius = 6371.0e3
    
    # Resolution
    nλ = 256
    nθ = 128
    nd = 1
    num_fourier = floor(Int64, nθ*(2/3))
    num_spherical = num_fourier + 1
    
    # Mesh
    mesh = Spectral_Spherical_Mesh(num_fourier = num_fourier, 
                                   num_spherical = num_spherical,
                                   nλ = nλ, 
                                   nθ = nθ, 
                                   nd = nd,
                                   radius = radius)
    
    # Coordinate
    cosθ = mesh.cosθ
    sinθ = mesh.sinθ
    
    # Main field
    grid_u = zeros(Float64, nλ, nθ, nd)
    grid_v = zeros(Float64, nλ, nθ, nd)
    grid_vor = zeros(Float64, nλ, nθ, nd)
    grid_div = zeros(Float64, nλ, nθ, nd)
    
    #
    if IC == "Case 1"
        # simple pure zonal-mean zonal flow, homogeneous
        for i = 1:nλ
            for j = 1:nθ
                for k = 1:nd
                    grid_u[i,j,k] = 10.0
                    grid_v[i,j,k] = 0.0
                    grid_vor[i,j,k] = 10.0*sinθ[j]/cosθ[j]/radius
                    grid_div[i,j,k] = 0.0
                end
            end
        end
    elseif IC == "Case 2"
        # simple pure zonal-mean zonal flow, centre on equator
        for i = 1:nλ
            for j = 1:nθ
                for k = 1:nd
                    grid_u[i,j,k] = 10.0*cosθ[j]^2
                    grid_v[i,j,k] = 0.0
                    grid_vor[i,j,k] = 10.0*3*sinθ[j]*cosθ[j]/radius
                    grid_div[i,j,k] = 0.0
                end
            end
        end
    elseif IC == "Case 3"
        # simple pure zonal-mean zonal flow, centre on poles, symmetric
        for i = 1:nλ
            for j = 1:nθ
                for k = 1:nd
                    grid_u[i,j,k] = 10.0*sinθ[j]^2
                    grid_v[i,j,k] = 0.0
                    grid_vor[i,j,k] = (10.0/radius/cosθ[j])*(2*sinθ[j]*cosθ[j]^2-sinθ[j]^3)
                    grid_div[i,j,k] = 0.0
                end
            end
        end
    elseif IC == "Case 4"
        # simple pure zonal-mean zonal flow, centre on poles, asymmetric
        for i = 1:nλ
            for j = 1:nθ
                for k = 1:nd
                    grid_u[i,j,k] = 10.0*sinθ[j]^3
                    grid_v[i,j,k] = 0.0
                    grid_vor[i,j,k] = 10.0*(3*sinθ[j]^2*cosθ[j]^2-sinθ[j]^4)/radius/cosθ[j]
                    grid_div[i,j,k] = 0.0
                end
            end
        end
    elseif IC == "Case 5"
        # simple pure zonal-mean meridional flow, homogeneous
        for i = 1:nλ
            for j = 1:nθ
                for k = 1:nd
                    grid_u[i,j,k] = 0.0
                    grid_v[i,j,k] = 10.0
                    grid_vor[i,j,k] = 0.0
                    grid_div[i,j,k] = -10.0*sinθ[j]/cosθ[j]/radius
                end
            end
        end
    elseif IC == "Case 6"
        # non-divergent flow, with multiple vortex
        λc = mesh.λc
        sinλ = sin.(λc)
        cosλ = cos.(λc)
        for i = 1:nλ
            for j = 1:nθ
                for k = 1:nd
                    grid_u[i,j,k] =  10.0*cosλ[i]*sinθ[j]
                    grid_v[i,j,k] = -10.0*sinλ[i]
                    grid_vor[i,j,k] = -20.0*cosλ[i]*cosθ[j]/radius
                    grid_div[i,j,k] = 0.0
                end
            end
        end
    elseif IC == "None"
        # non-divergent flow, with multiple vortex
        λc = mesh.λc
        θc = mesh.θc
        sinλ = sin.(λc)
        cosλ = cos.(λc)
        for i = 1:nλ
            for j = 1:nθ
                for k = 1:nd
                    grid_u[i,j,k] = 0.0
                    grid_v[i,j,k] = 0.0
                    grid_vor[i,j,k] = 4e-5*exp(-((θc[j]-19*pi/180)/(3*pi/180))^2)
                    grid_div[i,j,k] = 0.0
                end
            end
        end
        spe_vor = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd)
        spe_div = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd)
        Trans_Grid_To_Spherical!(mesh, grid_vor, spe_vor)
        Trans_Grid_To_Spherical!(mesh, grid_div, spe_div)
        UV_Grid_From_Vor_Div!(mesh, spe_vor, spe_div, grid_u, grid_v)
    else
        error("Initial condition ", IC, " has not implemented yet")
    end
    return mesh, grid_u, grid_v, grid_vor, grid_div
end


function UnitTest_CoordTrans(;IC::String = "Case 1",
                             iter_max::Int64 = 10)
    """
    TODO:: Add plots
    Test Trans_Grid_To_Spherical!, Trans_Spherical_To_Grid!
    """
    #
    mesh, grid_u, grid_v, grid_vor, grid_div = IC_Profiler(IC)
    nλ, nθ, nd, num_fourier, num_spherical = mesh.nλ, mesh.nθ, mesh.nd, mesh.num_fourier, mesh.num_spherical
    #
    grid_u_tmp = zeros(Float64, nλ, nθ, nd)
    grid_v_tmp = zeros(Float64, nλ, nθ, nd)
    grid_vor_tmp = zeros(Float64, nλ, nθ, nd)
    grid_div_tmp = zeros(Float64, nλ, nθ, nd)
    #
    grid_data = zeros(Float64, nλ, nθ, nd)
    grid_tmp = zeros(Float64, nλ, nθ, nd)
    spe_tmp = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd)
    
    # U
    grid_data .= grid_u
    grid_tmp .= grid_data
    for i = 1:iter_max
        Trans_Grid_To_Spherical!(mesh, grid_tmp, spe_tmp)
        Trans_Spherical_To_Grid!(mesh, spe_tmp, grid_tmp)
    end
    grid_u_tmp .= grid_tmp
    
    # V
    grid_data .= grid_v
    grid_tmp .= grid_data
    for i = 1:iter_max
        Trans_Grid_To_Spherical!(mesh, grid_tmp, spe_tmp)
        Trans_Spherical_To_Grid!(mesh, spe_tmp, grid_tmp)
    end
    grid_v_tmp .= grid_tmp
    
    # Vor
    grid_data .= grid_vor
    grid_tmp .= grid_data
    for i = 1:iter_max
        Trans_Grid_To_Spherical!(mesh, grid_tmp, spe_tmp)
        Trans_Spherical_To_Grid!(mesh, spe_tmp, grid_tmp)
    end
    grid_vor_tmp .= grid_tmp
    
    # Div
    grid_data .= grid_div
    grid_tmp .= grid_data
    for i = 1:iter_max
        Trans_Grid_To_Spherical!(mesh, grid_tmp, spe_tmp)
        Trans_Spherical_To_Grid!(mesh, spe_tmp, grid_tmp)
    end
    grid_div_tmp .= grid_tmp
    
    return (grid_u, grid_v, grid_vor, grid_div), 
           (grid_u_tmp, grid_v_tmp, grid_vor_tmp, grid_div_tmp), 
           mesh
end


function UnitTest_VelocityTrans(;IC::String = "Case 1",
                                iter_max::Int64 = 10)
    """
    TODO
    test different IC when given spectral?
    """
    mesh, grid_u, grid_v, grid_vor, grid_div  = IC_Profiler(IC)
    nλ, nθ, nd, num_fourier, num_spherical = mesh.nλ, mesh.nθ, mesh.nd, mesh.num_fourier, mesh.num_spherical
    
    # Vor & Div L2 error
    grid_u_tmp = zeros(Float64, nλ, nθ, nd)
    grid_v_tmp = zeros(Float64, nλ, nθ, nd)
    spe_vor_tmp = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd)
    spe_div_tmp = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd)
    
    grid_u_tmp .= grid_u
    grid_v_tmp .= grid_v
    
    for i = 1:iter_max
        Vor_Div_From_Grid_UV!(mesh, grid_u_tmp, grid_v_tmp, spe_vor_tmp, spe_div_tmp)
        UV_Grid_From_Vor_Div!(mesh, spe_vor_tmp, spe_div_tmp, grid_u_tmp, grid_v_tmp)
    end
    
    
    return (grid_u, grid_v), 
           (grid_u_tmp, grid_v_tmp), 
           mesh
end


function UnitTest_Advection(IC::String = "rigid_rotation")
    """
    """
    #
    mesh, grid_u, grid_v, grid_vor, grid_div  = IC_Profiler(IC)
    nλ, nθ, nd, num_fourier, num_spherical = mesh.nλ, mesh.nθ, mesh.nd, mesh.num_fourier, mesh.num_spherical 
    #
    cosθ, sinθ = mesh.cosθ, mesh.sinθ
    radius = mesh.radius
    λc = mesh.λc
    #
    grid_hs0 = zeros(Float64, nλ, nθ, nd)
    grid_∇hs0 = zeros(Float64, nλ, nθ, 2, nd)
    for k = 1:nd
        for i = 1:nλ
            for j = 1:nθ
                grid_hs0[i, j, k] = 25 * cosθ[j] - 30 * cosθ[j]^3 + 300 * sinθ[j]^2 * cosθ[j]^6 *sin(λc[i])
                grid_∇hs0[i, j, 1, k] = (300 * sinθ[j]^2 * cosθ[j]^5 *cos(λc[i]))/radius
                grid_∇hs0[i, j, 2, k] = (-25 * sinθ[j] + 90 * cosθ[j]^2*sinθ[j] + 600 * sinθ[j] * cosθ[j]^7*sin(λc[i]) - 1800 * sinθ[j]^3 * cosθ[j]^5*sin(λc[i]))/radius 
            end
        end
    end
    spe_hs0 = zeros(Complex{Float64}, num_fourier + 1, num_spherical + 1, nd)
    
    Trans_Grid_To_Spherical!(mesh, grid_hs0,  spe_hs0) 
    spe_cos_dλ_hs = zeros(Complex{Float64}, num_fourier + 1, num_spherical + 1, nd)
    spe_cos_dθ_hs = zeros(Complex{Float64}, num_fourier + 1, num_spherical + 1, nd)
    
    spe_dθ_hs0 = zeros(Complex{Float64}, num_fourier + 1, num_spherical + 1, nd)
    Trans_Grid_To_Spherical!(mesh, grid_∇hs0[:,:,2,:],  spe_dθ_hs0) 
            
    Compute_Gradient_Cos!(mesh, spe_hs0, spe_cos_dλ_hs, spe_cos_dθ_hs)
    
            
    spe_vor0 = zeros(Complex{Float64}, num_fourier + 1, num_spherical + 1, nd)
    spe_div0 = zeros(Complex{Float64}, num_fourier + 1, num_spherical + 1, nd)
    spe_vor1 = zeros(Complex{Float64}, num_fourier + 1, num_spherical + 1, nd)
    spe_div1 = zeros(Complex{Float64}, num_fourier + 1, num_spherical + 1, nd)
    grid_vor0 = zeros(Float64, nλ, nθ, nd)
    grid_div0 = zeros(Float64, nλ, nθ, nd)
    grid_vor1 = zeros(Float64, nλ, nθ, nd)
    grid_div1 = zeros(Float64, nλ, nθ, nd)
    
    grid_δhs = zeros(Float64, nλ, nθ, nd)
            
    Vor_Div_From_Grid_UV!(mesh, grid_u0, grid_v0, spe_vor0, spe_div0)
            
    Trans_Spherical_To_Grid!(mesh, spe_div0,  grid_div0)
    
    Vor_Div_From_Grid_UV!(mesh, grid_u0.*grid_hs0, grid_v0.*grid_hs0, spe_vor1, spe_div1)
    
    Trans_Spherical_To_Grid!(mesh, spe_div1,  grid_div1)
    
    Add_Horizontal_Advection!(mesh, spe_hs0, grid_u0, grid_v0, grid_δhs)
    grid_δhs_ref = -(grid_∇hs0[:, :, 1,:] .*grid_u0 + grid_∇hs0[:, :, 2,:] .* grid_v0)
    
    @show "adv approach error: ", norm(grid_δhs_ref - grid_δhs)
    @show "divergence approach error: ", norm(grid_δhs_ref + grid_div1 + grid_hs0.*grid_div0)
end


if abspath(PROGRAM_FILE) == @__FILE__
    # TODO: cosθ is poorly represented on spherical harmonic basis         
    test_derivative("rigid_rotation")
    test_derivative("random_1")
    test_derivative("random_2")
    test_derivative("mix")

    test_advection("rigid_rotation")
    test_advection("random_1")
    test_advection("random_2")
    test_advection("mix")
end