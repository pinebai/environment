#ifdef UNIFORM_PARTICLES
void assign_darkmatter_uniform(int last_dm_particle_id){
    int i, icell, level, ipart;
    double pos[nDim], vel[nDim];
    float Zsol = constants->Zsun;
    double pdt;
    int current_type;
    int current_id;
    
    current_id = last_dm_particle_id;
    current_type = 1;
    for(i=0;i<num_grid*num_grid*num_grid;i++)
	{
	    cell_center_position(i, pos);
	    icell = cell_find_position(pos);
	    vel[0] = 0;
	    vel[1] = 0;
	    vel[2] = 0;

	    level = cell_level(icell);
	    ipart = particle_alloc( current_id );
	    cart_assert( ipart > -1 && ipart < num_particles );
	    particle_x[ipart][0] = pos[0];
	    particle_x[ipart][1] = pos[1];
	    particle_x[ipart][2] = pos[2];
	    particle_v[ipart][0] = vel[0];
	    particle_v[ipart][1] = vel[1];
	    particle_v[ipart][2] = vel[2];

	    if( particle_species( current_id ) != current_type )
		{
		    cart_error("Assertion failed: particle_species(%d)=%d, current_type=%d",current_id,particle_species(current_id),current_type);
		}
	    particle_t[ipart] = tl[min_level];
	    particle_dt[ipart] = 0.0; /* dtl[level]; */
	    particle_mass[ipart] = particle_species_mass[particle_species(current_id)];

	    /* insert particle into cell linked list */
	    insert_particle( icell, ipart );
	    if(particle_x[ipart][0]>=num_grid || particle_x[ipart][0]<=0 ||
	       particle_x[ipart][1]>=num_grid || particle_x[ipart][1]<=0 ||
	       particle_x[ipart][2]>=num_grid || particle_x[ipart][2]<=0 ){
		cart_debug("%e %d %e ", particle_x[ipart][0],num_grid, particle_x[ipart][1]);
	    }

	    cart_assert(particle_x[ipart][0]<num_grid && particle_x[ipart][0]>0 );
	    cart_assert(particle_x[ipart][1]<num_grid && particle_x[ipart][1]>0 );
	    cart_assert(particle_x[ipart][2]<num_grid && particle_x[ipart][2]>0 );
	    
	    current_id++;
	}       
}
#endif /* UNIFORM_PARTICLES */

#ifdef UNIFORM_PARTICLES
    num_particle_species++; 
    i++;
    particle_species_mass[i] = uniform_particles_mass ;
    particle_species_indices[i+1] = particle_species_indices[i] + num_grid*num_grid*num_grid ; /* first index of specie */
    particle_species_num[i] = particle_species_indices[i+1] - particle_species_indices[i];
#endif
#ifdef UNIFORM_PARTICLES
double uniform_particles_mass;
#endif

////////////////
    if(iter_number == 0){
	cart_debug("assigning darkmatter");
	assign_darkmatter_model(fd, nhalo);

#ifdef UNIFORM_PARTICLES
	cart_assert(uniform_particles_mass>0);
	assign_darkmatter_uniform(nhalo);
#endif
//////////////


#ifdef UNIFORM_PARTICLES
    /* describe expected refinement */
    double code_tot_mass,dm1_mass,gas_mass,star_mass;
/*     const double model_gas_particle_mass =  8.96862e-05; /\* units 1e9*Msun*\/ */
/*     const double model_star_particle_mass =  4.484309e-04; */
/*     const double model_dm_particle_mass =  1.309238e-02; */
    /* the above should come from the read of ngas, nstars, nhalo*/

    code_tot_mass =  pow(num_grid,3.0);
    gas_mass  = ngas   * model_gas_particle_mass  * 1e9*constants->Msun/units->mass;
    dm1_mass  = nhalo   * model_dm_particle_mass * 1e9*constants->Msun/units->mass;
    star_mass = (nstars) * model_star_particle_mass * 1e9*constants->Msun/units->mass;

    uniform_particles_mass = (code_tot_mass - dm1_mass - gas_mass - star_mass)/pow(num_grid,3.0);
    
    double split_on_8gas,split_on_1dm;   
    split_on_8gas = 8 * model_gas_particle_mass*1e9*constants->Msun/units->mass;
    split_on_1dm = model_dm_particle_mass*1e9*constants->Msun/units->mass;
    cart_debug( "gas refinement=8x%e=%e ;  DM refinement=%e ; particles mass to place in every root cell =%e", 
		split_on_8gas/8.0, split_on_8gas,
		split_on_1dm,
		uniform_particles_mass
	);
#endif
