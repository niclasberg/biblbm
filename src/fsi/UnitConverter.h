#ifndef UNITCONVERTER_H_
#define UNITCONVERTER_H_

namespace plb {

namespace fsi {

template<class T>
class UnitConverter {
public:
	UnitConverter() : dx(1), dt(1), drho(1) { }

	/*
	 * Input:
	 *   dx = L_physical / L_lbm
	 *   dt = T_physical / T_lbm
	 *   drho = rho_physical / rho_lbm
	 */
	UnitConverter(T dx_, T dt_, T drho_)
	{
		init(dx_, dt_, drho_);
	}

	void set_si_to_lb_length_ratio(T dx_) { dx = dx_; }
	void set_si_to_lb_time_ratio(T dt_) { dt = dt_; }
	void set_si_to_lb_density_ratio(T drho_) { drho = drho_; }

	// Scalings
	T get_time_scale() const { return dt; }
	T get_length_scale() const { return dx; }
	T get_area_scale() const { return dx*dx; }
	T get_volume_scale() const { return dx*dx*dx;}

	T get_velocity_scale() const { return dx / dt; }
	T get_acceleration_scale() const { return dx / (dt*dt); }

	T get_density_scale() const { return drho; }
	T get_mass_scale() const { return drho * get_volume_scale(); }

	T get_force_scale() const { return get_mass_scale() * get_acceleration_scale(); }
	T get_force_density_scale() const { return get_density_scale() * get_acceleration_scale(); }
	T get_pressure_scale() const { return get_force_scale() / get_area_scale(); }

	T get_energy_scale() const { return get_force_scale() * get_length_scale(); }

	// lb to si unit conversions
	T lb_time_to_si(T time_lb) const { return time_lb * get_time_scale(); }
	T lb_length_to_si(T lb_length) const { return lb_length * get_length_scale(); }
	T lb_force_to_si(T lb_force) const { return lb_force * get_force_scale(); }
	T lb_force_density_to_si(T lb_force_density) const { return lb_force_density * get_force_density_scale(); }
	T lb_density_to_si(T dens_lb) const { return dens_lb * get_density_scale(); }
	T lb_velocity_to_si(T vel_lb) const { return vel_lb * get_velocity_scale(); }
	T lb_stress_to_si(T stress_lb) const { return stress_lb * get_force_scale() / get_area_scale(); }

	// si to lb unit conversion
	T si_time_to_lb(T time_si) const { return time_si / get_time_scale(); }
	T si_length_to_lb(T si_length) const { return si_length / get_length_scale(); }
	T si_force_to_lb(T si_force) const { return si_force / get_force_scale(); }
	T si_force_density_to_lb(T si_force_density) const { return si_force_density / get_force_density_scale(); }
	T si_pressure_to_lb(T si_pressure) const { return si_pressure / (get_force_scale() / get_area_scale()); }
	T si_density_to_lb(T dens_si) const { return dens_si / get_density_scale(); }
	T si_velocity_to_lb(T vel_si) const { return vel_si / get_velocity_scale(); }
	T si_stress_to_si(T stress_si) const { return stress_si * get_area_scale() / get_force_scale(); }
	T si_energy_to_lb(T energy_si) const { return energy_si / get_energy_scale(); }
	T si_force_per_length_to_lb(T force_per_length_si) const { return force_per_length_si * get_length_scale() / get_force_scale(); }

private:
	// Set properties
	T dt, dx, drho;
};

}

}




#endif /* UNITCONVERTER_H_ */
