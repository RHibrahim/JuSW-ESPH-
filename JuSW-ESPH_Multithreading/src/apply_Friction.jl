function apply_Friction!(int_fluid, vx, vy, h, hvx, hvy, fr_manning, d_hvx_dt, d_hvy_dt, grav)
    
    @inbounds @threads for i in int_fluid
        
        vx_magnitude = sqrt(vx[i]*vx[i] + vy[i]*vy[i])
        h_eff  = max(h[i], 1E-6)
		
		term = fr_manning[i]^2 * vx_magnitude / (h_eff^(4.0/3.0))
		
        d_hvx_dt[i] -= grav * hvx[i] * term
        d_hvy_dt[i] -= grav * hvy[i] * term

    end
    
end
