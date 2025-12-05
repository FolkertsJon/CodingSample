// This library is used to find the fast fourier transform (fft) of a data set from
// our oscilloscope data. Data sets going in need to have the structure of:
// time, ch1, ch2, ...
// where each entry above is a 1d vector<double> of the same length


#include<fftData.hh>

namespace User{
  // takes a data set with time, ch1, ch2, ... and returns frequency, fft1, fft2, ...
  std::vector<std::vector<double>> fftData(std::vector<std::vector<double>> dataSet){
    
    int nColumns = dataSet.size(); // get number of data columns
    int N = dataSet[0].size(); // get length of a data set
    double timeStep = dataSet[0][1]-dataSet[0][0];
    std::vector<std::vector<double>> fftData;

    if (nColumns < 2 ){
      std::cout << "This function requires the time base and at least one channel of data. There is only zero or one data set in what you've provided.\n";
      std::vector<double> negativeOne = {-1};
      fftData.push_back(negativeOne);
      return fftData;
    }


    // Calculate the frequency vector
    std::vector<double> freq(N/2+1);
    freq[0] = 0.0;
    for (int i = 1; i < N/2+1; ++i) {
      freq[i] = i/(timeStep*N);
    }
    fftData.push_back(freq);
    
    for( int i = 1; i < nColumns; i++){
      // Allocate memory for FFTW3 input and output arrays
      double *in = (double*) fftw_malloc(sizeof(double) * N);
      fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));
    
      // Create FFTW3 plan and execute the FFT
      fftw_plan plan = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
      for (int j = 0; j < N; j++) {
	in[j] = dataSet[i][j];
      }
      fftw_execute(plan);

      // Calculate the absolute magnitudes
      std::vector<double> magnitudes(N/2+1);
      //magnitudes[0] = std::abs(out[0][0])/N;
      for (int j = 0; j < N/2+1; j++) {
	magnitudes[j] = std::sqrt(std::pow(out[j][0]/N, 2) + std::pow(out[j][1]/N, 2));
      }
      fftData.push_back(magnitudes);
      fftw_destroy_plan(plan);
      fftw_free(in);
      fftw_free(out);
    }
    
    return fftData;
  }
}

