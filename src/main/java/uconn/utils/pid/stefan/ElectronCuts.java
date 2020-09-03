package uconn.utils.pid.stefan;

public class ElectronCuts {

    /**
     * Cherenkov nphe cut
     * @param nphe number of photoelectrons
     */
    public static boolean CC_nphe_cut(double nphe) {
        double nphe_min = 2;
        return nphe > nphe_min;
    }

    /**
     * PCAL energy
     * @param pcal_energy energy deposited in PCAL
     */
    public static boolean EC_outer_vs_EC_inner_cut(double pcal_energy) {
        final double edep_tight = 0.06, edep_medium = 0.07, edep_loose = 0.09;
        return pcal_energy > edep_medium;
    }


    /**
     * @param partp particle momentum
     * @param pcal_sector sector of hits in PCAL
     * @param pcal_energy energy in PCAL
     * @param ecin_energy energy in ECIN
     * @param ecout_energy energy in ECOUT
     */
    public static boolean EC_sampling_fraction_cut(double partp, int pcal_sector, double pcal_energy, double ecin_energy, double ecout_energy) {

        double[][] ecal_e_sampl_mu = {{  0.2531,  0.2550,  0.2514,  0.2494,  0.2528,  0.2521 },
            { -0.6502, -0.7472, -0.7674, -0.4913, -0.3988, -0.703  },
            {  4.939,  5.350,  5.102,  6.440,  6.149,  4.957  }
        };

        double[][] ecal_e_sampl_sigm = {{  2.726e-3,  4.157e-3,  5.222e-3,  5.398e-3,  8.453e-3,  6.533e-3 },
            {  1.062,  0.859,  0.5564,  0.6576,  0.3242,  0.4423   },
            { -4.089, -3.318, -2.078, -2.565, -0.8223, -1.274    }
        };

        double sigma_range = 3.5;

        double ectotal_energy = pcal_energy + ecin_energy + ecout_energy;
        int isec = pcal_sector-1;
        double mean = ecal_e_sampl_mu[0][isec] + ecal_e_sampl_mu[1][isec]/1000*Math.pow(partp-ecal_e_sampl_mu[2][isec],2);
        double sigma = ecal_e_sampl_sigm[0][isec] + ecal_e_sampl_sigm[1][isec]/(10*(partp-ecal_e_sampl_sigm[2][isec]));
        double upper_lim_total = mean + sigma_range * sigma;
        double lower_lim_total = mean - sigma_range * sigma;

        boolean pass_band = ectotal_energy/partp <= upper_lim_total && ectotal_energy/partp >= lower_lim_total;
        boolean pass_triangle = false;

        if(partp < 4.5) {
            pass_triangle = true;
        }
        else {
            pass_triangle = ecin_energy/partp > (0.2 -  pcal_energy/partp);
        }

        return pass_band && pass_triangle;
    }




    /**
     * @param pcal_sector sector of hits in PCAL
     * @param v pcal_lv
     * @param w pcal_lw
     */
    public static boolean EC_hit_position_fiducial_cut_homogeneous(int pcal_sector, double lv, double lw) {

        // Cut using the natural directions of the scintillator bars/ fibers:

        ///////////////////////////////////////////////////////////////////
        /// inbending:
        //
        double[] min_v_tight_inb = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
        double[] min_v_med_inb   = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
        double[] min_v_loose_inb = {9.0,  9.0,  9.0,  9.0,  9.0,  9.0 };
        //
        double[] max_v_tight_inb = {400, 400, 400, 400, 400, 400};
        double[] max_v_med_inb   = {400, 400, 400, 400, 400, 400};
        double[] max_v_loose_inb = {400, 400, 400, 400, 400, 400};
        //
        double[] min_w_tight_inb = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
        double[] min_w_med_inb   = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
        double[] min_w_loose_inb = {9.0,  9.0,  9.0,  9.0,  9.0,  9.0 };
        //
        double[] max_w_tight_inb = {400, 400, 400, 400, 400, 400};
        double[] max_w_med_inb   = {400, 400, 400, 400, 400, 400};
        double[] max_w_loose_inb = {400, 400, 400, 400, 400, 400};


        //////////////////////////////////////////////////////////////

        int isec = pcal_sector - 1;
        double min_v = min_v_loose_inb[isec];
        double max_v = max_v_loose_inb[isec];
        double min_w = min_w_loose_inb[isec];
        double max_w = max_w_loose_inb[isec];

        return lv > min_v && lv < max_v && lw > min_w && lw < max_w;
    }


    /**
     * @param dc_sector sector of hits in DC
     * @param region specify fiducial cuts for which region to use
     * @param x x for region 1 or 2 or 3 from REC::Traj
     * @param y y for region 1 or 2 or 3 from REC::Traj
     * @param partpid pid assigned to particle candidate
     * @param isinbending True if magnetic field is inbending
     */
    public static boolean DC_fiducial_cut_XY(int dc_sector, int region, double x, double y, int partpid, boolean isinbending) {

// new cut parameters for the linear cut based on x and y coordinates (inbending field):
// replace it in the function: bool DC_fiducial_cut_XY(int j, int region)
// (optimized for electrons, do not use it for hadrons)
//


        double[][][][] maxparams_in = {
            {   {{-14.563, 0.60032},{-19.6768, 0.58729},{-22.2531, 0.544896}},
                {{-12.7486, 0.587631},{-18.8093, 0.571584},{-19.077, 0.519895}},
                {{-11.3481, 0.536385},{-18.8912, 0.58099},{-18.8584, 0.515956}},
                {{-10.7248, 0.52678},{-18.2058, 0.559429},{-22.0058, 0.53808}},
                {{-16.9644, 0.688637},{-17.1012, 0.543961},{-21.3974, 0.495489}},
                {{-13.4454, 0.594051},{-19.4173, 0.58875},{-22.8771, 0.558029}}
            },
            {   {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
                {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}
            },
            {   {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
                {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}
            },
            {   {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
                {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}
            },
            {   {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
                {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}
            },
            {   {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
                {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}
            }
        };

        double[][][][] minparams_in = {
            {   {{12.2692, -0.583057},{17.6233, -0.605722},{19.7018, -0.518429}},
                {{12.1191, -0.582662},{16.8692, -0.56719},{20.9153, -0.534871}},
                {{11.4562, -0.53549},{19.3201, -0.590815},{20.1025, -0.511234}},
                {{13.202, -0.563346},{20.3542, -0.575843},{23.6495, -0.54525}},
                {{12.0907, -0.547413},{17.1319, -0.537551},{17.861, -0.493782}},
                {{13.2856, -0.594915},{18.5707, -0.597428},{21.6804, -0.552287}}
            },
            {   {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
                {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}
            },
            {   {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
                {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}
            },
            {   {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
                {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}
            },
            {   {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
                {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}
            },
            {   {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},
                {{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}},{{0, 0},{0, 0},{0, 0}}
            }
        };


        //double maxparams_out[6][6][3][2] =
        double[][][][] maxparams_out = {
        };


        //double minparams_out[6][6][3][2] =
        double[][][][] minparams_out = {
        };


        double[][][][] minparams = isinbending ? minparams_in : minparams_out;
        double[][][][] maxparams = isinbending ? maxparams_in : maxparams_out;

        double X = x;
        double Y = y;

        double X_new = X * Math.cos(Math.toRadians(-60*(dc_sector-1))) - Y * Math.sin(Math.toRadians(-60*(dc_sector-1)));
        Y = X * Math.sin(Math.toRadians(-60*(dc_sector-1))) + Y * Math.cos(Math.toRadians(-60*(dc_sector-1)));
        X = X_new;

        int pid = 0;

        switch (partpid)
        {
        case 11:
            pid = 0;
            break;
        case 2212:
            pid = 1;
            break;
        case 211:
            pid = 2;
            break;
        case -211:
            pid = 3;
            break;
        case 321:
            pid = 4;
            break;
        case -321:
            pid = 5;
            break;
        default:
            return false;
        }

        //if(inbending == true) pid = 0; // use only for electrons in inbending case

        double calc_min = minparams[pid][dc_sector - 1][region-1][0] + minparams[pid][dc_sector - 1][region-1][1] * X;
        double calc_max = maxparams[pid][dc_sector - 1][region-1][0] + maxparams[pid][dc_sector - 1][region-1][1] * X;

        return (Y > calc_min) && (Y < calc_max);
    }


    /**
     * @param pcal_sector sector of hits in PCAL
     * @param partvz Z in vertex
     */
    public static boolean DC_z_vertex_cut(int pcal_sector, double partvz, boolean isinbending) {

        double vz_min_sect_inb[] = {-13, -13, -13, -13, -13, -13};
        double vz_max_sect_inb[] = { 12,  12,  12,  12,  12,  12};

        double vz_min_sect_outb[] = {-18, -18, -18, -18, -18, -18};
        double vz_max_sect_outb[] = { 10,  10,  10,  10,  10,  10};

        double[] vz_min_sect = new double[6];
        double[] vz_max_sect = new double[6];

        for(int i = 0; i < 6; i++) {
            if(isinbending) {
                vz_min_sect[i] = vz_min_sect_inb[i];
                vz_max_sect[i] = vz_max_sect_inb[i];
            } else {
                vz_min_sect[i] = vz_min_sect_outb[i];
                vz_max_sect[i] = vz_max_sect_outb[i];
            }
        }

        int isec = pcal_sector - 1;
        double vz_min = vz_min_sect[isec];
        double vz_max = vz_max_sect[isec];

        return partvz > vz_min && partvz < vz_max;
    }

}

