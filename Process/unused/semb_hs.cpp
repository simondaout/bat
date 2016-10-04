/* Aufruf von Python :

    k  = Cm.otest (maxp, nostat, nsamp, ntimes, nstep, dimX, dimY, Gmint, new_frequence,
                   minSampleCount, latv, lonv, traveltime, traces)

  backSemb = otest (ncpus,nostat,nsamp,ntimes,nstep, dimX,dimY, mint, new_frequence, minSampleCount, latv, lonv, traveltimes, traces)

     Required arguments:
       ncpus  : input int
       nostat : input int
       nsamp  : input int
       ntimes : input int
       nstep  : input int
       dimX   : input int
       dimY   : input int
       mint   : input float
       new_frequence : input int
       minSampleCount : input int

       latv        : input rank-1 array('d') with bounds (dimX*dimY)
       lonv        : input rank-1 array('d') with bounds (dimX*dimY)
       traveltimes : input rank-1 array('d') with bounds (nostat*dimX*dimY)
       traces      : input rank-1 array('d') with bounds (nostat*minSampleCount)

     Return objects:
        backSemb : rank-1 array('d') with bounds (ntimes*dimX*dimY)
*/
/*  
  omp_set_num_threads(_sct);
  #pragma omp parallel shared (traveltime,trace,nstep,nostat,nsamp,ntimes,sembmax,latv,lonv,sembmaxX,sembmaxY) 
                       private(i,j,k,l,x,y) reduction(+:semb,nomin,denom,sum,relstart_samples)   
*/
def otest (ncpus, nostat, nsamp, ntimes, nstep, dimX,dimY, mint, new_frequence, minSampleCount, latv, lonv, traveltimes, traces) :
 
	for i in range (ntimes) :
      	
	      # +++++++++display data in files ++++ for each timestep one file +++++folders name same as event time++++
	      # timestep = boost::lexical_cast <double> (_origin->time().value().toString("%s"))+i*nstep; 
	      # cout << "TIMESTEP " << timestep <<endl;

	      sembpath = path+"/"+boost::lexical_cast<std::string>(i)+".ASC";
	      # cout << "SEMB " << sembpath << endl;

	      #sembout.open (sembpath.c_str(), ios::out);
	
	   # // if (!sembout) {	assert (false);}
	   # ++++++++++++++++++++++++write calculated data into file with header+++++++++++++++++++++++++++++++++++++++++++

		#sembout << "# "<<timestep <<" , "<<recordstarttime <<endl<< "# CalcWindowStart "
              #        <<_calcTimeWindow.startTime().toString("%F %T")<<endl<< "# CalcWindowEnd   "
              #        <<_calcTimeWindow.endTime().toString("%F %T")<<endl<< "# step: "<<_step<<"s| ntimes: "
              #        <<ntimes << "| winlen: "<<_winlen<<"s"<<endl;

		#sembout << "# ";

#//		for (unsigned int u = 0; u < sembstreams.size(); u++)  	{sembout <<sembstreams[u]<< " ";}

		#sembout << endl<< "# southwestlat: "<<Latul<<" dlat: "<<_gridspacing<<
		#	              " nlat: "<<_dimX<<endl<< "# southwestlon: "<<Lonul<<" dlon: "<<_gridspacing<<
		#		       " nlat: "<<_dimY<<endl<< "#ddepth: "<<"0"<<" ndepth: "<<"1"<<endl;

		#  loop over grid points
		sembmax  = 0;
		sembmaxX = 0;
		sembmaxY = 0;

		for j in range (dimX * dimY) :
		   semb  = 0;
	 	   nomin = denom = 0;
				
		   for l in range (nsamp) :
		 	sum = 0;	

			for k in range (nostat) :
                         relstart_samples = int ((traveltime[k][j] - mint) * new_frequence + 0.5) + i * nstep
			    sum   += trace[k][relstart_samples+l]
		           denom += trace[k][relstart_samples+l] * trace[k][relstart_samples+l]         			 
			}//for nostat

	 		nomin += sum * sum;
		   # endfor nsamp
	 		
	 	   x    = latv[j];
		   y    = lonv[j];
		   semb = nomin / (float (nostat) * denom);

		   #sembout <<x <<" "<< y<<" "<< semb << endl;
				
		   if semb > sembmax :
			sembmax  = semb;// search for maximum and position of maximum on semblance grid for given time step	   
			sembmaxX = latv[j];
			sembmaxY = lonv[j];				
		   #endif

	       #endfor dimX *dimY

	       # sembout << "# maximum semblance: "<<sembmax<<" at latitude: "
              #        << sembmaxX<<" at longitude: "<< sembmaxY<<endl;
	       # sembout.close();
		
	} // for ntimes
	
  } // pragma end	
