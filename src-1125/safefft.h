#pragma once


void initFFTWLock();
void freeFFTWLock();
void fft2d_safe(float* buf, int nsam);
void ifft2d_safe(float* buf, int nsam);
