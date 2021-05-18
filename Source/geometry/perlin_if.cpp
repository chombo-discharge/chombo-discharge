/*!
  @file perlin_if.cpp
  @brief Implementation of perlin_if.H
  @author Robert Marskar
  @date Sep. 2017
*/

#include "perlin_if.H"

#include "CD_NamespaceHeader.H"

//
int perlin_if::s_perlin[256] = { 151,160,137,91,90,15,
				 131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
				 190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
				 88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
				 77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
				 102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
				 135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
				 5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
				 223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
				 129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
				 251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
				 49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
				 138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180};

//
perlin_if::perlin_if(const Real     a_noiseAmp,
		     const RealVect a_noiseFreq,
		     const Real     a_persistence,
		     const int      a_octaves,
		     const bool     a_reseed){
#if CH_SPACEDIM < 2 || CH_SPACEDIM > 3
  MayDay::Abort("perlin_if::perlin_if - This will only compile in 2D or 3D");
#endif

  CH_assert(a_octaves >= 1);

  //
  m_noiseAmp    = a_noiseAmp;
  m_noiseFreq   = a_noiseFreq;
  m_persistence = a_persistence;
  m_octaves     = a_octaves;

  // Use Ken Perlin's original permutation table
  for (int i = 0; i < 256; i++){
    p[i]       = s_perlin[i];
    p[i + 256] = s_perlin[i];
  }


  // Reseed the permutation table
  if(a_reseed){
    reseed(); 
  }
}




//
perlin_if::perlin_if(const perlin_if& a_inputIF){

  //
  m_noiseAmp    = a_inputIF.m_noiseAmp;
  m_noiseFreq   = a_inputIF.m_noiseFreq;
  m_persistence = a_inputIF.m_persistence;
  m_octaves     = a_inputIF.m_octaves;
  for (int i = 0; i < 256; i++){
    p[i]           = a_inputIF.p[i];
    p[i + 256]     = a_inputIF.p[i + 256];
  }
}

//
perlin_if::~perlin_if(){
}

//
Real perlin_if::value(const RealVect& a_pos) const{
  return octaveNoise(a_pos);
}

//
BaseIF* perlin_if::newImplicitFunction() const{
  return static_cast<BaseIF*> (new perlin_if(*this));
}

//
int perlin_if::random (const int i) {
  return std::rand()%i;
}

//
void perlin_if::reseed(){

  //
  for (int i = 0; i < 256; i++){
    p[i] = i;
  }

  // Reseed the RNG
  //
  int seed = time(NULL);

  // Haven't had trouble with this, but I'm putting this here for safety. 
#ifdef CH_MPI
  int result = MPI_Bcast(&seed, 1, MPI_INT, 0, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){
    MayDay::Error("perlin_if::reseed() - broadcast failed");
  }
#endif
  std::srand(seed);

  
  
  random_shuffle(&p[0], &p[256], random);

  //
  for (int i = 0; i < 256; i++){
    p[i + 256] = p[i];
  }
}

double perlin_if::noise(const double a_x, const double a_y, const double a_z) const {

  // Lower cube corner
  const int X = (int) std::floor(a_x) & 255;
  const int Y = (int) std::floor(a_y) & 255;
  const int Z = (int) std::floor(a_z) & 255;

  // Relative distance wrt lower cube corner
  const double x = a_x - std::floor(a_x);
  const double y = a_y - std::floor(a_y);
  const double z = a_z - std::floor(a_z);

  // Fade curves
  const double u = fade(x);
  const double v = fade(y);
  const double w = fade(z);

  // Hash coordinates of 8 cube corners
  const int A  = p[X  ] + Y;
  const int AA = p[A  ] + Z;
  const int AB = p[A+1] + Z;
  const int B  = p[X+1] + Y;
  const int BA = p[B  ] + Z;
  const int BB = p[B+1] + Z;

  // Add blended results from 8 corners of cube
  return lerp(w, lerp(v, lerp(u, grad(p[AA  ], x  , y  , z   ),  // AND ADD
			      grad(p[BA  ], x-1, y  , z   )), // BLENDED
		      lerp(u, grad(p[AB  ], x  , y-1, z   ),  // RESULTS
			   grad(p[BB  ], x-1, y-1, z   ))),// FROM  8
	      lerp(v, lerp(u, grad(p[AA+1], x  , y  , z-1 ),  // CORNERS
			   grad(p[BA+1], x-1, y  , z-1 )), // OF CUBE
		   lerp(u, grad(p[AB+1], x  , y-1, z-1 ),
			grad(p[BB+1], x-1, y-1, z-1 ))));
}

Real perlin_if::noise(const RealVect& a_pos) const {

  Real x, y, z;
  x = a_pos[0];
  y = a_pos[1];
  z = 0.;
  
#if CH_SPACEDIM == 3
  z = a_pos[2];
#endif

  // Make noise [0, 1]
  return 0.5*noise(x, y, z) + 0.5;
}

Real perlin_if::octaveNoise(const RealVect& a_pos) const {
  Real result = 0.0;

  //
  Real normFrac = 0.;
  RealVect freq = m_noiseFreq;
  double amp    = 1.;//m_noiseAmp;

  // Add noise octaves
  for (int i = 0; i < m_octaves; ++i){
    result += noise(a_pos*freq)*amp;
    
    normFrac += amp;
    freq *= 1./m_persistence;
    amp  *= m_persistence;
  }

  // Normalize
  result *= m_noiseAmp;///normFrac;
  
  return result;
}

Real perlin_if::lerp(const Real t, const Real a, const Real b) const {
  return a + t * (b - a);
}

Real perlin_if::fade(const Real t) const {
  return t * t * t * (t * (t * 6 - 15) + 10); 
}

Real perlin_if::grad(const int hash, const double x, const double y, const double z) const {
  const int h    = hash & 15;
  const double u = h < 8 ? x : y;
  const double v = h < 4 ? y : h == 12 || h == 14 ? x : y;
  return ((h&1) == 0 ? u : -u) + ((h&2) == 0 ? v : -v);
}
#include "CD_NamespaceFooter.H"
