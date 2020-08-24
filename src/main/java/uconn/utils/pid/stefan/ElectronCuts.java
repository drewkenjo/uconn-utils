package uconn.utils.pid.stefan;

public class ElectronCuts {

    /**
     * PID cut
     * @param pid PDG code
     */
    public static boolean PID_cut(int pid) {
        return pid==11;
    }


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
     * @param region
     * @param x x for region 1 or 2 or 3 from REC::Traj
     * @param y y for region 1 or 2 or 3 from REC::Traj
     * @param partpid pid assigned to particle candidate
     * @param isinbending True if magnetic field is inbending
     */
    public static boolean DC_fiducial_cut_XY(int dc_sector, int region, double x, double y, int partpid, boolean isinbending) {

        //fitted values
        //double maxparams_in[6][6][3][2] =
        double[][][][] maxparams_in = {
            {   {{-7.49907, 0.583375},{-18.8174, 0.599219},{-23.9353, 0.574699}},
                {{-14.0547, 0.631533},{-14.4223, 0.597079},{-14.838, 0.547436}},
                {{-7.72508, 0.578501},{-18.7928, 0.56725},{-29.9003, 0.612354}},
                {{-6.12844, 0.566777},{-13.6772, 0.573262},{-26.1895, 0.591816}},
                {{-20.0718, 0.670941},{-9.4775, 0.511748},{-28.0869, 0.590488}},
                {{-9.52924, 0.591687},{-17.8564, 0.596417},{-23.5661, 0.576317}}
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

        //double minparams_in[6][6][3][2] =
        double[][][][] minparams_in = {
            {   {{7.62814, -0.56319},{18.2833, -0.587275},{20.2027, -0.54605}},
                {{9.20907, -0.586977},{10.493, -0.544243},{23.0759, -0.581959}},
                {{12.5459, -0.631322},{20.5635, -0.618555},{26.3621, -0.576806}},
                {{8.36343, -0.552394},{14.7596, -0.554798},{29.5554, -0.60545}},
                {{16.3732, -0.663303},{10.0255, -0.533019},{31.6086, -0.617053}},
                {{8.20222, -0.567211},{20.0181, -0.605458},{22.2098, -0.567599}}
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
            {   {{-15.0486, 0.72351},{-19.2202, 0.611078},{-27.4513, 0.535668}},
                {{-27.5984, 0.856463},{-16.7873, 0.601081},{-23.1163, 0.552932}},
                {{-10.4537, 0.698218},{-13.8137, 0.559658},{-35.1087, 0.605458}},
                {{-15.2615, 0.747703},{-18.0476, 0.598169},{-26.6229, 0.540454}},
                {{-13.2849, 0.669102},{-11.4523, 0.50548},{-14.619, 0.448236}},
                {{-13.9748, 0.706217},{-19.2459, 0.61829},{-32.0969, 0.586403}}
            },
            {   {{1.01396, 0.498685},{-15.5912, 0.576023},{-15.1113, 0.532929}},
                {{-5.78163, 0.558017},{-16.2599, 0.579884},{-18.757, 0.553068}},
                {{-5.85198, 0.563657},{-12.8175, 0.560827},{-16.0893, 0.53454}},
                {{-4.73676, 0.550713},{-15.9589, 0.578277},{-16.4839, 0.533793}},
                {{-6.246, 0.567365},{-14.435, 0.570721},{-15.4366, 0.531068}},
                {{-6.97001, 0.576377},{-14.6786, 0.568779},{-19.7683, 0.560218}}
            },
            {   {{-5.11874, 0.542303},{-15.0805, 0.574983},{-14.4774, 0.525724}},
                {{-12.0383, 0.614845},{-13.7267, 0.563175},{-15.0412, 0.526371}},
                {{-6.03209, 0.557031},{-13.7801, 0.565532},{-17.3829, 0.544538}},
                {{-5.11456, 0.544146},{-12.2437, 0.550863},{-17.309, 0.543019}},
                {{-10.4872, 0.598711},{-10.4236, 0.532708},{-14.3734, 0.53108}},
                {{-1.0749, 0.520489},{-12.8321, 0.562083},{-16.5056, 0.5455}}
            },
            {   {{-4.02259, 0.524186},{-17.1253, 0.594224},{-18.2492, 0.548083}},
                {{-9.8687, 0.599888},{-16.4612, 0.586687},{-13.8826, 0.523052}},
                {{-6.09038, 0.565961},{-17.7112, 0.599618},{-15.4107, 0.527006}},
                {{-7.39969, 0.589083},{-17.9769, 0.601895},{-11.6827, 0.508671}},
                {{-9.4989, 0.597539},{-14.5596, 0.570631},{-11.0431, 0.511108}},
                {{-4.78634, 0.541011},{-15.7703, 0.58688},{-13.4413, 0.527737}}
            },
            {   {{-8.52794, 0.597096},{-15.7673, 0.578371},{-15.3708, 0.532371}},
                {{-9.93026, 0.599944},{-18.9386, 0.600393},{-16.9471, 0.538744}},
                {{-6.97351, 0.580954},{-14.787, 0.571151},{-14.7902, 0.535692}},
                {{-10.1381, 0.61534},{-15.7761, 0.578828},{-12.4601, 0.516824}},
                {{-6.96808, 0.575769},{-16.8735, 0.59029},{-16.5667, 0.546479}},
                {{-4.66291, 0.557931},{-12.8325, 0.564619},{-15.4321, 0.532078}}
            },
            {   {{-7.04376, 0.566423},{-17.1298, 0.584139},{-5.45416, 0.378063}},
                {{-8.44714, 0.587986},{-18.1342, 0.594626},{1.15458, 0.235862}},
                {{-5.94346, 0.565299},{-13.1923, 0.544233},{-9.35298, 0.356599}},
                {{-6.21599, 0.565299},{-17.3319, 0.592404},{-9.93697, 0.490374}},
                {{-4.59187, 0.547947},{-18.6261, 0.617462},{-4.56516, 0.318076}},
                {{-4.51481, 0.54217},{-15.5234, 0.586091},{-3.33085, 0.406969}}
            }
        };


        //double minparams_out[6][6][3][2] =
        double[][][][] minparams_out = {
            {   {{10.4534, -0.646633},{20.3057, -0.634014},{19.8023, -0.456707}},
                {{8.39762, -0.604501},{16.3439, -0.594578},{24.3968, -0.550637}},
                {{11.7765, -0.682773},{18.4467, -0.616169},{33.5937, -0.587674}},
                {{10.3321, -0.651872},{20.2277, -0.636832},{18.6805, -0.5109}},
                {{7.23971, -0.624146},{15.3243, -0.578226},{19.042, -0.50646}},
                {{11.3041, -0.656507},{17.9646, -0.598628},{20.7513, -0.490774}}
            },
            {   {{2.38471, -0.519764},{17.37, -0.58695},{18.0912, -0.543571}},
                {{5.54502, -0.55162},{15.7411, -0.574163},{14.124, -0.527875}},
                {{-2.07845, -0.477191},{17.6555, -0.581425},{19.773, -0.545526}},
                {{0.0286526, -0.510017},{14.3589, -0.569907},{18.4092, -0.553169}},
                {{4.65436, -0.546023},{16.5496, -0.579403},{16.1278, -0.53216}},
                {{2.42162, -0.520676},{13.5759, -0.557783},{15.4258, -0.523538}}
            },
            {   {{5.42714, -0.547252},{17.9145, -0.586462},{12.5151, -0.50246}},
                {{1.34431, -0.50385},{15.4816, -0.56814},{14.8332, -0.530269}},
                {{4.21559, -0.534567},{16.107, -0.56405},{13.833, -0.514593}},
                {{4.20743, -0.536347},{15.2913, -0.571454},{16.8639, -0.534427}},
                {{6.1354, -0.558451},{15.8315, -0.575873},{16.3529, -0.541569}},
                {{8.68088, -0.575933},{16.4584, -0.576239},{15.6839, -0.534505}}
            },
            {   {{4.68623, -0.547042},{15.677, -0.565814},{15.9174, -0.525793}},
                {{5.63431, -0.550919},{16.3397, -0.582825},{13.0206, -0.519814}},
                {{10.7146, -0.609206},{17.798, -0.588528},{17.5752, -0.535447}},
                {{5.60915, -0.55906},{15.9808, -0.582329},{15.4179, -0.531727}},
                {{5.02979, -0.550918},{16.5882, -0.588454},{14.8774, -0.521648}},
                {{11.0076, -0.616694},{16.3595, -0.578062},{14.1671, -0.52539}}
            },
            {   {{7.23685, -0.571469},{12.413, -0.523034},{10.2567, -0.43924}},
                {{2.89157, -0.524202},{13.7971, -0.541647},{16.9268, -0.51881}},
                {{7.39136, -0.564481},{15.7659, -0.579543},{13.2347, -0.51334}},
                {{5.30074, -0.555581},{8.29014, -0.455451},{10.9363, -0.504937}},
                {{8.93705, -0.586499},{15.6892, -0.575975},{10.5961, -0.455508}},
                {{9.26221, -0.591506},{14.466, -0.565013},{15.2475, -0.53195}}
            },
            {   {{7.31173, -0.576571},{15.7123, -0.594739},{19.0354, -0.499379}},
                {{8.07107, -0.591361},{17.1351, -0.603906},{-0.239706, -0.24005}},
                {{7.83093, -0.577219},{11.579, -0.516213},{14.2709, -0.394561}},
                {{9.06302, -0.600711},{18.5039, -0.616063},{-3.86975, -0.213802}},
                {{4.71687, -0.547772},{18.0103, -0.604119},{10.7901, -0.437388}},
                {{7.17946, -0.57369},{11.2732, -0.397794},{-2.79229, -0.137468}}
            }
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

