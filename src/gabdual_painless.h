#ifndef GABDUAL_PAINLESS_H__
#define GABDUAL_PAINLESS_H__

/** Compute the first dl samples of the Gabor frame operator diagonal
 *
 * \param[in]   g    Original window
 * \param[in]  gl    Length of the windows
 * \param[in]   a    Hop factor
 * \param[in]   M    Number of channels
 * \param[in]  dl    Number of frame diagonal samples
 * \param[out]  d    Frame diagonal
 */
void gabframediag(const double *g, int gl, int a, int M, int dl, double *d);

/** Compute canonical dual window for painless Gabor system
 *
 * \param[in]   g    Original window
 * \param[in]  gl    Length of the windows
 * \param[in]   a    Hop factor
 * \param[in]   M    Number of channels
 * \param[out] gd    Canonical dual window
 */
void gabdual_painless(const double g[], int gl, int a, int M, double gd[]);

#endif  // GABDUAL_PAINLESS_H__