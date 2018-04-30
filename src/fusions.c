#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "iio.h"
#include "pickopt.h"
#include "fusions.h"


//double weighted_mean_with_binary_weight(
//        double *values,
//        double nv,
//        double *weight,
//        double weight_power,
//        double *binary_weight
//        )
//{
//    double weighted_sum = 0;
//    double weights_sum = 0;
//    for (int i = 0; i < nv; i++)
//    {
//        weighted_sum += !isnan(values[i]) * !isnan(binary_weight[i]) * pow(weight[i],
//                weight_power) * values[i];
//        weights_sum += !isnan(values[i]) * !isnan(binary_weight[i]) * pow(weight[i],
//                weight_power);
//    }
//    return (weighted_sum / weights_sum);
//}
//
//
//double weighted_mean(
//        double *values,
//        double nv,
//        double *weight,
//        double weight_power
//        )
//{
//    double weighted_sum = 0;
//    double weights_sum = 0;
//    for (int i = 0; i < nv; i++)
//    {
//        weighted_sum += !isnan(values[i]) * pow(weight[i], weight_power) *
//            values[i];
//        weights_sum += !isnan(values[i]) * pow(weight[i],
//                weight_power);
//    }
//    return (weighted_sum / weights_sum);
//}
//
//double weighted_mean_with_binary_weight_for_median(
//        double *values,
//        double nv,
//        double *weight,
//        double weight_power,
//        double *binary_weight,
//        double tmp_med
//        )
//{
//    double weighted_sum = 0;
//    double weights_sum = 0;
//    for (int i = 0; i < nv; i++)
//    {
//        weighted_sum += !isnan(values[i]) * !isnan(binary_weight[i]) * pow(weight[i],
//                weight_power) * values[i] / fabs(values[i]-tmp_med);
//        weights_sum += !isnan(values[i]) * !isnan(binary_weight[i]) * pow(weight[i],
//                weight_power) / fabs(values[i] - tmp_med);
//    }
//    return (weighted_sum / weights_sum);
//}
//
//double weighted_mean_for_median(
//        double *values,
//        double nv,
//        double *weight,
//        double weight_power,
//        double tmp_med
//        )
//{
//    double weighted_sum = 0;
//    double weights_sum = 0;
//    for (int i = 0; i < nv; i++)
//    {
//        weighted_sum += !isnan(values[i]) * pow(weight[i], weight_power) *
//            values[i] / fabs(values[i]-tmp_med);
//        weights_sum += !isnan(values[i]) * pow(weight[i], weight_power) /
//            fabs(values[i] - tmp_med);
//    }
//    return (weighted_sum / weights_sum);
//}
//
//double weighted_median_with_binary_weight(
//        double *values,
//        double nv,
//        double *weight,
//        double weight_power,
//        double *binary_weight
//        )
//{
//    double median = 0;
//    for (int i = 0; i < 3; i++)           // converge en peu d'itérations
//        median = weighted_mean_with_binary_weight_for_median(values, nv, weight, weight_power, binary_weight, median); 
//    return median;
//}
//
//double weighted_median(
//        double *values,
//        double nv,
//        double *weight,
//        double weight_power
//        )
//{
//    double median = 0;
//    for (int i = 0; i < 3; i++)           // converge en peu d'itérations
//        median = weighted_mean_for_median(values, nv, weight, weight_power, median); 
//    return median;
//}
//
//double p_error_to_minimise(
//        double *values,
//        double nv,
//        double m,
//        double p
//        )
//{
//    double sum = 0;
//    for (int i = 0; i < nv; i++)
//        sum += pow(fabs(values[i] - m), p); 
//    return pow(sum/nv, 1/p);
//}
//
//double frechet_p_centroid(
//        double *values,
//        double nv,
//        double p)
//{
//    double centroid = 0;
//    for (int i = 0; i < nv; i++)
//        if (p_error_to_minimise(values, nv, centroid, p) 
//                > p_error_to_minimise(values, nv, values[i], p))
//            centroid = values[i];
//    return centroid;
//}
//
//double weighted_p_error_to_minimise(
//        double *values,
//        double nv,
//        double *weight,
//        double weight_power,
//        double m,
//        double p
//        )
//{
//    double sum = 0;
//    double weights_sum = 0;
//    for (int i = 0; i < nv; i++){
//        sum += !isnan(values[i]) * pow(weight[i], weight_power) *
//            pow(fabs(values[i] - m), p); 
//        weights_sum += !isnan(values[i]) * pow(weight[i], weight_power);}
//    return pow(sum/weights_sum, 1/p);
//}
//
//double weighted_frechet_p_centroid(
//        double *values,
//        double nv,
//        double *weight,
//        double weight_power,
//        double p)
//{
//    double centroid = 0;
//    for (int i = 0; i < nv; i++)
//        if (weighted_p_error_to_minimise(values, nv, weight, weight_power,
//                    centroid, p) 
//                > weighted_p_error_to_minimise(values, nv, weight,  weight_power,
//                    values[i], p))
//            centroid = values[i];
//    return centroid;
//}
//
//double weighted_p_error_to_minimise_with_binary_weights(
//        double *values,
//        double nv,
//        double *weight,
//        double weight_power,
//        double *binary_weight,
//        double m,
//        double p
//        )
//{
//    double sum = 0;
//    double weights_sum = 0;
//    for (int i = 0; i < nv; i++){
//        sum += !isnan(values[i]) * !isnan(binary_weight[i]) * pow(weight[i],
//                weight_power) * pow(fabs(values[i] - m), p);
//        weights_sum += !isnan(values[i]) * !isnan(binary_weight[i]) * pow(weight[i],
//                weight_power);}
//    return pow(sum/weights_sum, 1/p);
//}
//
//double weighted_frechet_p_centroid_with_binary_weights(
//        double *values,
//        double nv,
//        double *weight,
//        double weight_power,
//        double *binary_weight,
//        double p)
//{
//    double centroid = 0;
//    for (int i = 0; i < nv; i++)
//        if (weighted_p_error_to_minimise_with_binary_weights(values, nv, weight,
//                    weight_power, binary_weight, centroid, p) >
//                weighted_p_error_to_minimise_with_binary_weights(values, nv, weight,
//                    weight_power, binary_weight, values[i], p))
//            centroid = values[i];
//    return centroid;
//}
//
//double fusion_for_one_vertex(
//        double *values,
//        double *weight,
//        double weight_power,
//        double *binary_weight,
//        int method,
//        double p,
//        int h)
//{
//    switch (method){
//        case 10:
//            return weighted_mean(values, h, weight, weight_power);
//        case 11:
//            return weighted_mean_with_binary_weight(values, h, weight,
//                    weight_power, binary_weight);
//        case 20:
//            return weighted_median(values, h, weight, weight_power);
//        case 21:
//            return weighted_median_with_binary_weight(values, h, weight,
//                    weight_power, binary_weight);
//        case 30:
//            return weighted_frechet_p_centroid(values, h, weight, weight_power,
//                    p);
//        case 31:
//            return weighted_frechet_p_centroid_with_binary_weights(values, h,
//                    weight, weight_power, binary_weight, p);
//        default : 
//            fprintf(stderr, "Invalid fusion method\n");
//            return NAN;
//    }
//}
//
//void fusion_for_entire_vectors(
//        double *rgb,
//        double *scalars,
//        double weight_power,
//        double *sun,
//        int method,
//        double frechet_p,
//        int  w, int h,
//        double *out)
//{
//    double values[h];
//    double weight[h];
//    double binary_weight[h];
//    for (int i = 0; i < w; i++)
//        for (int k = 0; k < 3; k++)
//        {
//            for (int j = 0; j < h; j++)
//            {
//                values[j] = rgb[3 * (i + j*w) + k];
//                weight[j] = scalars[3 * (i + j*w) + k];
//                binary_weight[j] = sun[3 * (i + j*w) + k];
//            }
//            out[3*i+k] = fusion_for_one_vertex(values, weight, weight_power, binary_weight, method, frechet_p, h);
//        }
//}


int main(int c, char *v[])
{
    double p = atof(pick_option(&c, &v, "p", "0"));
    if (c < 5)
        return fprintf(stderr, "usage:\n\t"
                "%s all_rgb.tif all_scalars.tif all_real_sun.tif all_theoretical_sun.tif"
                "weight_power, method, out.tif\n"
                "available methods :\n"
                " - 10 weighted mean\n"
                " - 11 weighted mean with binary weights\n"
                " - 20 weighted median\n"
                " - 21 weighted median with binary weights\n"
                " - 30 weighted frechet centroid\n"
                " - 31 weighted frechet centroid with binary weights\n"
                " if using frechet centroid add p variable -p : \n"
                , *v);

    char *filename_rgb = v[1];
    char *filename_scalars = v[2];       // scaled mesh
    char *filename_real_sun = v[3];
    char *filename_theoretical_sun = v[4];        // panchromatic satellite image .ntf
    double weight_power = atof(v[5]);
    int method  = atoi(v[6]);
    char *filename_out = v[7];

    // read vertices utm coordinates
    int w, h, pd;
//    double *rgb = iio_read_image_double_vec(filename_rgb, &w, &h, &pd);
//    if (!rgb)
//        return fprintf(stderr, "iio_read(%s) failed\n", filename_rgb);
//
//    // read vertices utm coordinates
//    double *real_sun = iio_read_image_double_vec(filename_real_sun, &w, &h, &pd);
//    if (!real_sun)
//        return fprintf(stderr, "iio_read(%s) failed\n", filename_real_sun);
//
//    // read vertices utm coordinates
//    double *theoretical_sun = iio_read_image_double_vec(filename_theoretical_sun, &w, &h, &pd);
//    if (!theoretical_sun)
//        return fprintf(stderr, "iio_read(%s) failed\n", filename_theoretical_sun);
//
//    // read vertices utm coordinates
//    double *scalars = iio_read_image_double_vec(filename_scalars, &w, &h, &pd);
//    if (!scalars)
//        return fprintf(stderr, "iio_read(%s) failed\n", filename_scalars);
//
//    double *out = malloc(3 * w * sizeof(double));
//
//    fusion_for_entire_vectors(rgb, scalars, weight_power, real_sun, method, p, w, h, out);
//
//    // save outputs
//    iio_save_image_double_vec(filename_out, out, w, 1, 3);

    // free allocated memory
//    free(out); 

    return 0;
}












































