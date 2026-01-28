function update_openboundary!(int_openboundary, wse, h, vx, vy, z_b, volume, itype, MAX_NB_FLUID, IDbc, numOpenBound, numData_OpenBound, timeOpenBound, h_timeOpenBound, vx_timeOpenBound, time, boundary_normals_X, boundary_normals_Y, Openboundary_neighbor, Openboundary_Wij)
    
    h_TempOpenBound  = zeros(Float64, numOpenBound)
    vx_TempOpenBound = zeros(Float64, numOpenBound)

    # Linear interpolation of h_TempOpenBound and vx_TempOpenBound at time
    @inbounds @threads for i in 1:numOpenBound  # loop over boundaries

        htmp  = 0.0
        vxtmp = 0.0

        ndata = numData_OpenBound[i]

        @inbounds for ii = 1:(ndata - 1)
            t1 = timeOpenBound[ii, i]
            t2 = timeOpenBound[ii + 1, i]

            if time <= t2
                α = (time - t1) / (t2 - t1)

                htmp  = h_timeOpenBound[ii, i]  + α * (h_timeOpenBound[ii + 1, i]  - h_timeOpenBound[ii, i])
                vxtmp = vx_timeOpenBound[ii, i] + α * (vx_timeOpenBound[ii + 1, i] - vx_timeOpenBound[ii, i])
                break
            end
        end

        h_TempOpenBound[i]  = htmp
        vx_TempOpenBound[i] = vxtmp
        
    end
    
    @inbounds @threads for i in int_openboundary

        local_vx_OpenBound1 = 0.0
        local_vx_OpenBound2 = 0.0
        local_h_OpenBound   = 0.0
        local_z_b_OpenBound = 0.0
        local_sum_wij       = 0.0

        @inbounds for k = 1:MAX_NB_FLUID
            
            j = Openboundary_neighbor[k, i]
            j == 0 && continue
            
            w    = Openboundary_Wij[k, i]

            coef = volume[j] * w

            # Accumulate contributions from fluid particles (j) to boundary particle (i)
            local_vx_OpenBound1 += coef * vx[j]
            local_vx_OpenBound2 += coef * vy[j]
            local_h_OpenBound   += coef * h[j]
            local_z_b_OpenBound += coef * z_b[i]
            local_sum_wij       += coef
            
        end

        # Loop to normalize and apply boundary conditions 
        if local_sum_wij > 0.0 

            vx_OpenBound1 = local_vx_OpenBound1 / local_sum_wij
            vx_OpenBound2 = local_vx_OpenBound2 / local_sum_wij
            h_OpenBound   = local_h_OpenBound   / local_sum_wij
            z_b_OpenBound = local_z_b_OpenBound / local_sum_wij

        else

            vx_OpenBound1 = 0.0
            vx_OpenBound2 = 0.0
            h_OpenBound   = h_TempOpenBound[IDbc[i]]  # Default to boundary condition values
            z_b_OpenBound = 0.0

        end

        # Transform velocities to the new coordinate system
        upi_p = + vx_OpenBound1 * boundary_normals_X[i] + vx_OpenBound2 * boundary_normals_Y[i]  # Normal velocity
        vpi_p = - vx_OpenBound1 * boundary_normals_Y[i] + vx_OpenBound2 * boundary_normals_X[i]  # Tangential velocity

        if itype[i] == 2  # Inflow boundary
            
            up_ob_p = vx_TempOpenBound[IDbc[i]]     # Prescribed inflow velocity
            vp_ob_p = 0.0                           # Zero tangential velocity at inflow

            if up_ob_p > sqrt(grav*h_TempOpenBound[IDbc[i]])    # Supercritical flow
                h[i] = h_TempOpenBound[IDbc[i]]                 
            else                                     # Subcritical flow
                sqh = (0.5 / sqrt(grav)) * (up_ob_p - upi_p) + sqrt(h_OpenBound)
                h[i] = sqh * sqh
            end

        elseif itype[i] == -2 # Outflow boundary

            if upi_p > sqrt(grav*h_OpenBound)     # Supercritical flow
                h[i]    = h_OpenBound
                up_ob_p = upi_p
                vp_ob_p = vpi_p
            else                                     # Subcritical flow
                h[i]    = h_TempOpenBound[IDbc[i]]
                up_ob_p = upi_p + 2 * sqrt(grav) * (sqrt(h_OpenBound) - sqrt(h_TempOpenBound[IDbc[i]]))
                vp_ob_p = vpi_p
            end

        end

        # Transform velocities back to original coordinate system
        vx[i] = + up_ob_p * boundary_normals_X[i] - vp_ob_p * boundary_normals_Y[i]
        vy[i] = + up_ob_p * boundary_normals_Y[i] + vp_ob_p * boundary_normals_X[i]
   
        z_b[i] = z_b_OpenBound
        wse[i] = h[i] + z_b[i]
        
    end

    return nothing

end