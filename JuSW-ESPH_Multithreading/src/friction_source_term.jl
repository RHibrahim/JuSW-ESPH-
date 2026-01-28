function apply_Friction!(particle_ids, vx, h, hvx, fr_manning, d_hvx_dt)
    
    @inbounds @threads for i in particle_ids
        
        vx1 = vx[i, 1]
        vx2 = vx[i, 2]

        vx_magnitude = sqrt(vx1*vx1 + vx2*vx2)
        h_eff  = max(h[i], 1E-6)
		
		term = fr_manning[i]^2 * vx_magnitude / (h_eff^(4.0/3.0))
		
        d_hvx_dt[i, 1] -= grav * hvx[i, 1] * term
        d_hvx_dt[i, 2] -= grav * hvx[i, 2] * term

    end
    
end
