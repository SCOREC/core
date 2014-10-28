#ifndef MILO_H
#define MILO_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct milo* milo_t;

milo_t milo_new(int w, int h, const char* servname, int port);
void milo_free(milo_t m);

/* points are arrays of 3 doubles (xyz),
   colors are arrays of 3 doubles (rgb), where each color is from 0 to 1 */

void milo_clear(milo_t m, double* background);
void milo_dot(milo_t m, double* point, double* color);
void milo_line(milo_t m, double* point_a, double* point_b, double* color);
void milo_triangle(milo_t m, double* point_a, double* point_b, double* point_c,
    double* color);
void milo_text(milo_t m, double* point, const char* text, double* color);

void milo_zoom(milo_t m, double factor);
void milo_pan(milo_t m, double x, double y);
void milo_spin(milo_t m, double radians);
void milo_tilt(milo_t m, double radians);

void milo_run(milo_t m);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
