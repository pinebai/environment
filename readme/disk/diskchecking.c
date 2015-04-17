
// after remap_star_ids in the iter_numer ==0


////////////////////////////////////////////////////// from tree_debug.c
/* 	int count=0; */
/* 	int species_count_total[100]; */
/*         int species_count[100]; */
/* 	int i; */

/*         for ( i = 0; i < num_particle_species; i++ ) { */
/* 	    species_count[i] = 0; */
/* 	    species_count_total[i] = 0; */
/*         } */
/*         for ( i = 0; i < num_particles; i++ ) { */
/* 	    if ( particle_level[i] != FREE_PARTICLE_LEVEL ) { */
/* 		if( particle_id[i] >= num_particles_total ) */
/* 		    { */
/* 			cart_debug("Incorrect particle[%d] id=%d, num_particles_total=%d",i,particle_id[i],num_particles_total); */
/* 		    } */
/* 		count++; */
/* 		species_count[ particle_species( particle_id[i] ) ]++; //snl   */
/* 	    } */
/*         } */

/* 	MPI_Allreduce( species_count, species_count_total, num_particle_species, MPI_INT, MPI_SUM, mpi.comm.run ); */
/* 	MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_INT, MPI_SUM, mpi.comm.run ); */

/* 	cart_debug(" count %d particles total %d ",count,num_particles_total);  */
/*         for ( i = 0; i < num_particle_species; i++ ) { //snl */
/* 	    cart_debug( "spec %d counted %d supposedly there %d",i, species_count_total[i], particle_species_num[i] ); */
/* 	    cart_assert( species_count_total[i] == particle_species_num[i] ); */
/*         } */
//	exit(1);
//////////////////////////////////////////////////////
