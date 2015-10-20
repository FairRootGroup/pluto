double PData:: sampleM(const int &id) {
    if ( makeStaticData()->getParticleTotalWidth(id)<0.0003 || 
	 !getDepth(id,1) ) // do not sample narrow resonances
	return makeStaticData()->getParticleMass(id); 



