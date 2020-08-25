package uconn.utils.pid.stefan;

public class HadronCuts {

    /**
     * DC fiducial cut for hadrons
     * @param dc_sector sector of hits in DC
     * @param region specify fiducial cuts for which region to use
     * @param trajx x for region 1 or 2 or 3 from REC::Traj
     * @param trajy y for region 1 or 2 or 3 from REC::Traj
     * @param trajz z for region 1 or 2 or 3 from REC::Traj
     * @param partpid pid assigned to particle candidate
     * @param isinbending True if magnetic field is inbending
     */
    public static boolean DC_fiducial_cut_theta_phi(int dc_sector, int region, double trajx, double trajy, double trajz, int partpid, boolean isinbending) {

        //fitted values inbending:

        double[][][][] maxparams_in = {
            {   {{-35.1716, 25.102, -0.750281, 5.34679e-05},{-39.1633, 28.5551, -1.13429, 0.00419047},{-33.7705, 24.8068, -0.811239, 0.00138345}},
                {{-36.2389, 26.7979, -1.08147, 0.0050898},{-43.643, 31.6783, -1.49203, 0.00872922},{-54.4042, 40.6516, -2.52393, 0.0205649}},
                {{-38.3238, 26.1667, -0.777077, 0.000264835},{-34.2011, 24.2843, -0.696392, 3.75866e-12},{-36.4636, 25.8712, -0.786592, 2.24421e-10}},
                {{-31.8019, 23.154, -0.653992, 2.69968e-05},{-34.6637, 24.6043, -0.714901, 2.02675e-10},{-36.7209, 26.2469, -0.828638, 0.000340435}},
                {{-33.4016, 24.6901, -0.779889, 0.000430557},{-35.4583, 24.7491, -0.707953, 2.18559e-10},{-37.7335, 28.1547, -1.1986, 0.00582395}},
                {{-34.7808, 24.6988, -0.719936, 5.73299e-10},{-54.5797, 40.9138, -2.57493, 0.0213354},{-38.4972, 28.3142, -1.21741, 0.00640373}}
            },
            {   {{-2.25358e-08, 12.631, -0.767619, 0.00739811},{-8.09501, 15.9098, -0.844083, 0.00667995},{-1.48113e-06, 12.2061, -0.73167, 0.0074309}},
                {{-2.10872e-07, 12.6689, -0.765156, 0.00720044},{-4.88862, 14.0376, -0.687202, 0.00506307},{-4.59793e-06, 11.5553, -0.591469, 0.00536957}},
                {{-1.13504e-08, 12.6011, -0.746025, 0.00687498},{-6.97884, 15.1788, -0.765889, 0.00570532},{-1.29468, 12.3844, -0.667561, 0.00619226}},
                {{-2.91953e-09, 13.883, -0.999624, 0.0104257},{-4.9855, 13.8864, -0.661348, 0.0048371},{-2.29438e-08, 11.8341, -0.668486, 0.00669247}},
                {{-2.02824e-08, 13.3855, -0.91158, 0.00926769},{-3.29092e-08, 10.8294, -0.382323, 0.00178367},{-4.59027e-06, 11.9414, -0.663872, 0.00625769}},
                {{-3.73322e-09, 12.6126, -0.723548, 0.0062217},{-4.56248, 14.1574, -0.727805, 0.00560108},{-2.39381e-08, 12.0663, -0.6651, 0.00602544}}
            },
            {   {{-1.45923e-08, 13.0297, -0.828302, 0.00795271},{-5.41905, 13.2753, -0.503236, 0.00255607},{-3.67719, 12.1358, -0.462905, 0.00308219}},
                {{-9.953e-10, 11.549, -0.52816, 0.00378771},{-8.47154, 15.9863, -0.826166, 0.0062936},{-6.43715, 13.9081, -0.618535, 0.0046102}},
                {{-4.68458e-08, 12.9481, -0.781613, 0.00689754},{-3.46617, 12.2786, -0.440121, 0.00205448},{-4.43519, 10.9372, -0.210059, 3.69283e-10}},
                {{-4.18414e-07, 13.1542, -0.811251, 0.00714402},{-4.63166, 13.7769, -0.657207, 0.0047586},{-1.99278e-05, 11.3993, -0.575232, 0.00532141}},
                {{-7.07189e-10, 13.2814, -0.88476, 0.00874389},{-5.08373, 14.4384, -0.750795, 0.00586116},{-6.9642e-05, 9.50651, -0.189316, 3.07274e-06}},
                {{-5.85515e-08, 12.5116, -0.688741, 0.00557297},{-1.86306, 11.985, -0.482567, 0.00279836},{-4.94295e-07, 10.1342, -0.316715, 0.00176254}}
            },
            {   {{-0.0157256, 11.1508, -0.415185, 0.00186904},{-13.6561, 19.4418, -1.15773, 0.00989432},{-6.24969e-07, 10.5776, -0.329325, 0.00103488}},
                {{-2.5686e-08, 11.4797, -0.476772, 0.00264288},{-0.0475099, 10.1207, -0.244786, 3.13032e-06},{-4.6875e-07, 12.019, -0.63598, 0.00543214}},
                {{-0.00702545, 11.1294, -0.407207, 0.00171263},{-7.27687, 15.5, -0.807858, 0.0062086},{-5.15078, 12.6368, -0.348584, 9.2687e-12}},
                {{-8.14106e-08, 13.28, -0.818164, 0.00703758},{-7.60722, 14.4871, -0.588662, 0.00326244},{-1.70764e-06, 12.0413, -0.63961, 0.00541784}},
                {{-1.09281, 11.5573, -0.41311, 0.00155228},{-3.71599, 12.8335, -0.521472, 0.00296792},{-0.000410815, 12.4833, -0.72999, 0.0066601}},
                {{-0.652641, 12.2766, -0.554202, 0.00314615},{-8.42824, 15.5087, -0.710609, 0.00447051},{-14.9692, 21.5885, -1.47528, 0.0136615}}
            },
            {   {{-5.58945, 17.4004, -1.34516, 0.0142099},{-14.9585, 20.4538, -1.25118, 0.0106617},{-12.0069, 16.4545, -0.727162, 0.00495418}},
                {{-7.03048, 17.3519, -1.1831, 0.0111308},{-7.30641, 15.8503, -0.850952, 0.00648446},{-10.2549, 15.6139, -0.648352, 0.00380506}},
                {{-9.73111e-09, 13.498, -0.932479, 0.00939708},{-8.38053, 15.5588, -0.711323, 0.00433827},{-12.3097, 16.6403, -0.741362, 0.0050708}},
                {{-7.38905, 17.2652, -1.15517, 0.0109165},{-1.11835e-07, 10.4637, -0.301972, 0.000612754},{-12.2182, 17.4958, -0.919555, 0.00747512}},
                {{-0.492676, 14.4148, -1.0959, 0.0116708},{-5.34309, 14.3258, -0.691954, 0.00480109},{-12.5443, 16.1047, -0.59594, 0.00280171}},
                {{-4.08375e-07, 12.2846, -0.655278, 0.00525956},{-8.93101, 16.4266, -0.861853, 0.00644623},{-11.8406, 17.0417, -0.826301, 0.00596028}}
            },
            {   {{-9.29415, 16.5566, -0.831923, 0.00562504},{-0.954483, 10.5813, -0.265766, 3.24615e-05},{-6.87423, 14.892, -0.76495, 0.00639603}},
                {{-18.8913, 19.3123, -0.711917, 0.00227889},{-13.9788, 18.5678, -0.940183, 0.00664397},{-11.7696, 18.3415, -1.04368, 0.0083506}},
                {{-3.82873, 12.7727, -0.425968, 0.000789835},{-9.81221, 14.6531, -0.471092, 0.00131406},{-14.2392, 15.9895, -0.430525, 2.20712e-12}},
                {{-1.76975e-07, 11.4006, -0.420134, 0.00141302},{-3.11764, 10.9707, -0.245823, 2.23044e-12},{-17.6005, 22.2881, -1.39992, 0.0117791}},
                {{-0.767518, 11.6824, -0.456275, 0.00214005},{-5.28047, 12.65, -0.350658, 9.80081e-05},{-0.0888832, 11.508, -0.49197, 0.00301269}},
                {{-4.72388, 15.8507, -1.00574, 0.00876768},{-2.80649, 11.4056, -0.301812, 0.000190262},{-13.0484, 18.665, -1.08614, 0.00960977}}
            }
        };

        double[][][][] minparams_in = {
            {   {{37.289, -27.5201,1.12866, -0.00526111},{45.3103, -33.5226,1.72923, -0.0114495},{61.5709, -47.6158,3.4295, -0.0316429}},
                {{36.6259, -27.4064,1.16617, -0.00604629},{50.3751, -37.5848,2.19621, -0.0169241},{35.1563, -26.514,1.09795, -0.00545864}},
                {{27.2367, -20.3068,0.517752, -0.000335432},{39.0489, -28.6903,1.24306, -0.0065226},{41.0208, -30.0339,1.30776, -0.00626721}},
                {{29.261, -21.7041,0.613556, -0.000774652},{39.5304, -29.1388,1.34116, -0.00823818},{44.5313, -33.4056,1.77581, -0.0123965}},
                {{36.5659, -25.119,0.714074, -2.65397e-11},{31.6524, -22.6934,0.613977, -5.46634e-10},{34.7312, -24.9901,0.749061, -1.22922e-09}},
                {{33.154, -23.8803,0.685794, -1.13236e-10},{42.6731, -31.0799,1.40425, -0.00730816},{46.4732, -35.6988,2.10144, -0.0164771}}
            },
            {   {{2.40451, -15.0848,1.05504, -0.0103356},{8.93825, -16.5995,0.925874, -0.00767902},{7.23439e-08, -12.5963,0.814574, -0.00864749}},
                {{6.2953e-07, -12.6365,0.732206, -0.00639165},{12.6866, -18.7831,1.0952, -0.00923029},{3.12805e-07, -12.5395,0.795535, -0.00828991}},
                {{2.69495, -14.8778,1.00751, -0.00975373},{6.05446, -14.6778,0.767457, -0.00636729},{3.94741e-07, -11.1038,0.524109, -0.00471514}},
                {{2.31558e-07, -11.5073,0.494316, -0.00303611},{5.66995, -14.5948,0.740956, -0.00561851},{4.40475e-06, -9.57062,0.20354, -0.000213213}},
                {{2.74277e-08, -13.3573,0.886651, -0.00857992},{9.98849e-05, -11.524,0.531486, -0.00391441},{8.50811e-07, -9.72224,0.240264, -0.000781498}},
                {{6.9021e-08, -11.8859,0.53864, -0.00325092},{10.0169, -16.9153,0.921593, -0.00752414},{9.90518e-07, -11.9578,0.697029, -0.00717645}}
            },
            {   {{6.87966e-10, -12.8497,0.757379, -0.00651612},{16.2087, -19.3776,0.951508, -0.00645029},{14.513, -18.8625,1.05542, -0.00918985}},
                {{1.07197e-07, -12.5469,0.703086, -0.00585238},{0.0871522, -9.22628,0.159628, -0.000343326},{12.1181, -17.5575,0.940249, -0.00788125}},
                {{2.10191e-09, -12.2782,0.661926, -0.00555279},{12.5105, -17.9998,0.951807, -0.00732845},{12.8043, -17.8322,0.972401, -0.00841528}},
                {{8.11926e-10, -12.7225,0.737941, -0.00647355},{7.50649, -15.987,0.889398, -0.00729282},{0.174541, -10.0266,0.306882, -0.00186093}},
                {{3.81202e-09, -12.0926,0.598943, -0.00430458},{8.72368, -17.2511,1.06348, -0.00953327},{1.5205, -9.86713,0.183806, -6.40377e-12}},
                {{1.37378e-07, -12.9247,0.769722, -0.00664936},{8.53877, -16.6167,0.946138, -0.00788745},{8.47417, -14.3897,0.581492, -0.00387111}}
            },
            {   {{2.50079e-07, -12.5209,0.678491, -0.00528954},{12.6171, -18.4205,1.01802, -0.00807702},{10.4903, -18.0981,1.10546, -0.00971519}},
                {{5.87069e-07, -12.0075,0.585538, -0.00416654},{11.1348, -17.5468,0.943652, -0.00729083},{0.949201, -10.5869,0.267536, -6.04802e-05}},
                {{1.14857, -11.1478,0.345528, -0.000841836},{10.9482, -17.1647,0.909605, -0.00722404},{8.7569e-08, -10.4446,0.316302, -0.00101964}},
                {{1.09759e-06, -11.5019,0.48435, -0.00277852},{0.637937, -10.7065,0.316211, -0.000801127},{5.67144e-07, -12.88,0.831252, -0.00835441}},
                {{1.68853, -11.2582,0.308152, -7.81686e-12},{9.44238, -17.1892,1.00561, -0.00864837},{1.20713e-07, -12.2246,0.669321, -0.0057622}},
                {{0.00217558, -10.8858,0.347928, -0.000790679},{11.8583, -17.6423,0.923581, -0.00703041},{3.24078, -13.4024,0.668777, -0.00504175}}
            },
            {   {{6.04158, -16.8155,1.13335, -0.0105359},{8.24786, -17.0204,1.05097, -0.00941875},{11.7617, -17.202,0.864472, -0.00649032}},
                {{3.70947, -13.0663,0.513818, -0.00222627},{16.7022, -21.9618,1.42869, -0.012705},{6.8993, -14.8192,0.740813, -0.00585407}},
                {{2.18472e-06, -11.9461,0.583354, -0.00423414},{6.51489e-07, -10.5669,0.353028, -0.00166977},{12.5113, -16.5038,0.709888, -0.00471964}},
                {{0.812719, -11.3245,0.390183, -0.00134086},{2.97251, -11.9374,0.338592, -4.36096e-13},{13.8844, -17.5707,0.818446, -0.00581811}},
                {{1.55496, -14.4569,0.949497, -0.00857237},{0.34359, -10.5041,0.286497, -0.000346977},{14.4141, -18.7457,1.01652, -0.00845189}},
                {{1.26317e-08, -11.1424,0.434251, -0.00236267},{6.58119, -15.8546,0.930324, -0.00801288},{4.41865, -11.1991,0.234652, -7.43723e-10}}
            },
            {   {{6.87926, -12.8949,0.334733, -6.38494e-06},{35.2336, -32.2007,2.21489, -0.020555},{6.80949, -16.8945,1.19056, -0.0127558}},
                {{0.95782, -12.4625,0.599979, -0.00405342},{20.4051, -23.1936,1.42408, -0.0120792},{10.277, -16.1457,0.785186, -0.00612069}},
                {{0.236196, -11.6165,0.458613, -0.002018},{12.8771, -19.6785,1.26163, -0.0115917},{5.21194e-08, -12.551,0.78718, -0.00794713}},
                {{8.40778, -14.9001,0.534967, -0.00147246},{15.9376, -20.9945,1.2908, -0.0110556},{10.4773, -16.2238,0.783386, -0.00593478}},
                {{3.21187, -12.1221,0.348938, -8.70415e-14},{13.8983, -19.1128,1.04727, -0.00797426},{11.6342, -18.8428,1.18853, -0.0107619}},
                {{3.7311, -12.4292,0.419345, -0.00134704},{6.92884, -13.2494,0.391862, -0.000767396},{5.5939, -14.4175,0.729195, -0.00568477}}
            }
        };


        //fitted values outbending
        double[][][][] maxparams_out = {
            {   {{-3.69457, 12.3755, -0.41328, 0.00129631},{-54.3237, 40.3308, -2.39952, 0.0181339},{-39.8661, 27.1428, -0.907303, 0.00220974}},
                {{-37.6199, 26.2865, -0.826366, 0.000862203},{-72.4212, 54.7953, -4.04856, 0.0373308},{-21.1791, 17.0759, -0.391795, 0.00151085}},
                {{-0.421685, 10.482, -0.272111, 8.69408e-05},{-43.3635, 32.746, -1.6541, 0.0101454},{-62.6387, 41.1869, -1.97298, 0.0107022}},
                {{-42.0766, 29.6387, -0.993426, 1.97101e-09},{-44.7036, 33.0587, -1.64131, 0.0099416},{-47.2703, 32.6109, -1.46533, 0.00817871}},
                {{-22.2035, 20.6894, -0.689051, 0.000592423},{-74.6572, 54.7065, -3.83999, 0.0351952},{-38.9183, 25.7212, -0.711499, 2.5796e-12}},
                {{-52.078, 45.571, -3.71942, 0.0376577},{-65.4047, 49.1723, -3.36623, 0.0288435},{-53.9611, 35.9294, -1.58589, 0.00772417}}
            },
            {   {{-2.20312e-07, 13.0916, -0.864184, 0.0086342},{-6.44026e-08, 12.056, -0.675801, 0.00643464},{-20.2596, 23.5977, -1.545, 0.0141047}},
                {{-4.42537e-05, 10.2799, -0.322454, 0.00154825},{-1.63659e-07, 11.0228, -0.451412, 0.00308633},{-8.5382, 15.6903, -0.785315, 0.00602734}},
                {{-2.32088, 11.6343, -0.363509, 0.000902217},{-0.301128, 12.0319, -0.643794, 0.00581994},{-22.4378, 25.2772, -1.73656, 0.0164181}},
                {{-7.40627, 13.601, -0.382439, 2.45262e-05},{-5.50415e-08, 11.9792, -0.652368, 0.00597647},{-15.1608, 20.6455, -1.33827, 0.0127123}},
                {{-0.203913, 10.7032, -0.322123, 0.000691162},{-1.73184e-07, 10.735, -0.379993, 0.00196037},{-0.155443, 10.1794, -0.249841, 6.24278e-05}},
                {{-1.87352e-07, 12.4226, -0.730141, 0.0068049},{-1.40236e-07, 12.5356, -0.750615, 0.00719921},{-16.8681, 21.8555, -1.43078, 0.0131935}}
            },
            {   {{-8.89326e-08, 10.0681, -0.240869, 9.9612e-12},{-15.2705, 21.635, -1.55291, 0.0166645},{-10.5976, 17.9928, -1.08432, 0.00950807}},
                {{-0.00389562, 10.2092, -0.254082, 4.15737e-06},{-9.16032e-11, 10.527, -0.334641, 0.00129061},{-9.63013e-07, 11.0668, -0.42453, 0.0022955}},
                {{-2.40163e-06, 13.4151, -0.949883, 0.0107662},{-1.60937e-07, 10.5128, -0.35046, 0.00173787},{-29.2647, 30.1252, -2.20552, 0.0213809}},
                {{-2.69733e-08, 11.7703, -0.589854, 0.00482124},{-3.77564e-08, 11.3764, -0.527037, 0.00416671},{-4.85047, 13.7737, -0.650441, 0.0047428}},
                {{-3.90816e-07, 12.2683, -0.692591, 0.00625884},{-9.70203e-10, 11.0335, -0.438323, 0.00275342},{-2.54193, 13.5404, -0.76861, 0.00684486}},
                {{-3.23439e-10, 10.7412, -0.348557, 0.00113794},{-1.79623, 11.7499, -0.449432, 0.00247294},{-13.1393, 19.4689, -1.17148, 0.00984086}}
            },
            {   {{-5.07611e-08, 11.7796, -0.516966, 0.00295389},{-4.87018, 12.2727, -0.322719, 9.12315e-06},{-35.9369, 31.015, -1.95133, 0.0169834}},
                {{-1.32385e-07, 11.6454, -0.495467, 0.00272602},{-2.70664, 12.0151, -0.434014, 0.00203292},{-8.97137, 15.0453, -0.646138, 0.00429196}},
                {{-7.92247e-09, 12.5189, -0.682231, 0.00539531},{-0.0942499, 10.3465, -0.280521, 0.000405358},{-19.7485, 21.7919, -1.24334, 0.0105088}},
                {{-8.50093e-11, 10.739, -0.302295, 5.6862e-11},{-0.184771, 10.4358, -0.285869, 0.000389546},{-21.9469, 24.9675, -1.77893, 0.0183075}},
                {{-4.34589, 12.5902, -0.362849, 4.996e-15},{-0.000684493, 10.6055, -0.332363, 0.00104632},{-21.328, 22.0864, -1.20993, 0.00989151}},
                {{-0.0202168, 12.0097, -0.539165, 0.00299034},{-0.5239, 10.7167, -0.309141, 0.000535617},{-10.0299, 16.3179, -0.812315, 0.00617078}}
            },
            {   {{-0.169908, 10.902, -0.353938, 0.00100715},{-3.2818, 13.2193, -0.65495, 0.00515117},{-0.013532, 8.51331, -0.070239, 1.755e-05}},
                {{-8.51985e-08, 11.6512, -0.56808, 0.00453582},{-1.2381e-07, 10.6653, -0.368149, 0.00181989},{-9.30287e-08, 10.0352, -0.254321, 0.000417053}},
                {{-0.150407, 10.6338, -0.308676, 0.000481694},{-0.00186321, 10.4259, -0.303092, 0.00073092},{-21.3328, 28.0803, -2.37912, 0.025101}},
                {{-14.4411, 19.817, -1.13705, 0.00894685},{-6.25263e-09, 11.7414, -0.586098, 0.00478932},{-5.49193, 16.1248, -1.11306, 0.0115644}},
                {{-1.54761, 12.0015, -0.462506, 0.00204729},{-5.72883, 14.9638, -0.795325, 0.00616222},{-50.229, 45.8456, -3.88803, 0.0414729}},
                {{-40.7531, 33.6269, -2.03771, 0.01609},{-1.33363e-09, 11.9894, -0.614358, 0.004924},{-27.2506, 29.2602, -2.1426, 0.0203235}}
            },
            {   {{-1.62999e-10, 14.0422, -1.03609, 0.0107179},{-6.71565, 15.6964, -0.887791, 0.00740777},{-38.9148, 32.9935, -2.09023, 0.0177295}},
                {{-1.09078e-05, 13.4131, -0.878092, 0.00825152},{-15.0102, 21.6968, -1.4935, 0.0138851},{-19.5261, 20.3932, -0.969464, 0.00661531}},
                {{-1.39619e-08, 12.3593, -0.618488, 0.00415536},{-5.38271e-07, 11.5631, -0.512607, 0.00334452},{-23.0902, 24.7093, -1.57315, 0.0140132}},
                {{-1.73908e-08, 12.0348, -0.591608, 0.00423834},{-8.35134, 17.3066, -1.11555, 0.010407},{-2.74909e-07, 9.59202, -0.216455, 0.000527479}},
                {{-0.0449157, 10.5243, -0.334389, 0.00134555},{-0.0143489, 10.0993, -0.2434, 1.57595e-10},{-22.3661, 23.2499, -1.32946, 0.0108047}},
                {{-5.83731e-07, 14.5234, -1.14022, 0.0122177},{-1.4586e-08, 11.6946, -0.520935, 0.00324975},{-12.4252, 16.3216, -0.652566, 0.00365791}}
            }
        };

        double[][][][] minparams_out = {
            {   {{3.73672, -12.3584,0.390616, -0.000795415},{51.644, -37.8546,1.99228, -0.0119973},{32.3551, -22.9742,0.624096, -4.30811e-05}},
                {{6.11614, -13.6358,0.491668, -0.0018637},{47.5098, -35.902,1.97535, -0.0134876},{82.9536, -58.2741,4.12662, -0.0378612}},
                {{0.000950108, -7.99619,0.000506416, -0.0020788},{64.0688, -47.8642,3.16007, -0.025878},{70.0064, -50.3249,3.38975, -0.029639}},
                {{37.0145, -35.0316,2.61892, -0.0250306},{14.5954, -15.6554,0.426733, -0.000879865},{28.9035, -21.5279,0.610475, -0.00087271}},
                {{5.65685, -13.3347,0.400781, -1.46612e-11},{67.3504, -50.152,3.33677, -0.0270726},{47.0772, -32.1506,1.38851, -0.00719898}},
                {{8.95987, -15.1646,0.585477, -0.00246174},{41.6154, -29.7967,1.1817, -0.00403765},{61.1631, -41.6465,2.32522, -0.0175271}}
            },
            {   {{8.80954e-10, -11.0364,0.413853, -0.00210254},{6.50072e-08, -11.2505,0.501571, -0.00380973},{10.9643, -17.4701,0.989297, -0.00860789}},
                {{2.33292e-08, -11.2353,0.470728, -0.00309666},{2.29373e-07, -11.2458,0.50218, -0.00383969},{29.5429, -29.9965,2.19166, -0.021366}},
                {{1.61826e-08, -11.861,0.577321, -0.00433276},{2.9436e-07, -11.5738,0.581015, -0.00503307},{19.5142, -23.451,1.58724, -0.0151339}},
                {{2.07231e-09, -12.7453,0.751184, -0.00664181},{1.77802e-07, -11.4574,0.537367, -0.00422656},{12.5683, -18.4632,1.05475, -0.00892182}},
                {{7.6216e-08, -13.9769,1.01051, -0.0107372},{1.33092e-08, -11.9128,0.628521, -0.00550105},{13.5537, -20.1708,1.32578, -0.0123213}},
                {{9.25941, -19.658,1.51566, -0.0157124},{6.25983e-10, -11.6806,0.599263, -0.00532588},{17.0479, -22.0046,1.47474, -0.0140475}}
            },
            {   {{4.65436e-08, -11.1925,0.466196, -0.00308992},{18.4968, -22.5122,1.4594, -0.0135962},{18.9488, -23.3348,1.57414, -0.0146183}},
                {{3.67722e-08, -10.9985,0.428395, -0.00257574},{16.3745, -21.0105,1.3093, -0.0119156},{11.4404, -18.6679,1.15919, -0.010306}},
                {{1.46846e-08, -10.865,0.398638, -0.00212392},{20.7337, -23.3738,1.46852, -0.0130115},{28.2098, -28.9406,2.05908, -0.0197782}},
                {{0.237058, -10.4694,0.271985, -1.08731e-07},{2.32759, -11.9354,0.469887, -0.00291497},{13.287, -20.8621,1.49656, -0.0148999}},
                {{0.000149907, -10.4632,0.294713, -0.000431947},{6.96663, -15.3946,0.845078, -0.00724722},{11.0939, -17.4733,0.944239, -0.00747728}},
                {{3.10006e-08, -10.1416,0.247764, -1.36913e-11},{5.41915, -14.6085,0.795369, -0.00684375},{5.89127, -13.0881,0.453024, -0.0020325}}
            },
            {   {{4.16588e-09, -12.9305,0.749425, -0.00611725},{5.65263, -14.1661,0.637395, -0.00400239},{4.66325, -12.9519,0.565753, -0.00442033}},
                {{8.0428e-08, -13.1625,0.836744, -0.00778246},{12.3243, -18.8718,1.11103, -0.00917354},{7.20312, -16.0935,0.987223, -0.00930883}},
                {{0.00147165, -10.4992,0.280542, -1.79846e-06},{3.20232, -11.6892,0.350774, -0.00101099},{8.14117e-08, -10.9813,0.524839, -0.00507885}},
                {{0.470888, -13.5446,0.820782, -0.00768941},{3.9697, -13.0821,0.540847, -0.00303209},{3.44817, -12.3932,0.533804, -0.00414144}},
                {{1.05038e-08, -10.6539,0.297078, -6.04694e-05},{15.0983, -21.1791,1.38383, -0.0124058},{17.3666, -20.3986,1.16663, -0.0102393}},
                {{8.49365e-07, -13.765,0.964056, -0.00956575},{9.38084, -16.7385,0.904339, -0.00707907},{12.1048, -17.3704,0.91318, -0.00757461}}
            },
            {   {{10.6378, -19.5017,1.45275, -0.017057},{1.24368e-08, -10.5134,0.338985, -0.00143696},{37.3291, -35.1606,2.60092, -0.0242728}},
                {{19.1614, -24.0851,1.73932, -0.0185466},{14.1293, -19.8382,1.21613, -0.0107037},{20.9629, -24.0839,1.60283, -0.015173}},
                {{0.000450804, -8.15062,0.0103867, -2.00709e-05},{5.72496, -14.338,0.717819, -0.00567964},{16.9428, -21.8075,1.4216, -0.0131736}},
                {{6.15991e-10, -11.5278,0.536105, -0.00402223},{2.17842e-07, -10.5338,0.327427, -0.0010898},{20.7387, -24.3028,1.65004, -0.0155857}},
                {{0.650351, -10.6177,0.275393, -6.4664e-08},{8.05811, -16.1558,0.913735, -0.00788487},{0.308897, -10.2816,0.275186, -0.000561299}},
                {{0.427836, -10.168,0.240458, -5.90042e-06},{2.30661, -12.8686,0.664796, -0.00562626},{0.00499667, -11.6585,0.62597, -0.00619261}}
            },
            {   {{9.01249e-07, -11.8437,0.494125, -0.00223452},{14.3941, -21.2365,1.46048, -0.0137349},{13.7095, -15.4704,0.408961, -0.000312145}},
                {{0.000251044, -11.3084,0.438545, -0.0020791},{0.00847078, -12.6769,0.804431, -0.00836705},{1.09388, -9.66797,0.175278, -1.8721e-11}},
                {{4.04693e-10, -11.9001,0.585913, -0.00440376},{5.05178, -12.1514,0.31134, -0.000112735},{30.8105, -28.0795,1.73625, -0.0151639}},
                {{3.86607e-11, -13.471,0.889111, -0.0083617},{8.86591e-09, -9.25745,0.163052, -6.08491e-12},{27.1358, -24.3255,1.23326, -0.00891886}},
                {{0.196086, -11.7392,0.480055, -0.00224614},{0.18667, -10.5859,0.287231, -6.53153e-06},{14.8865, -17.1338,0.653576, -0.00333176}},
                {{2.7955e-07, -13.1311,0.848222, -0.00812719},{29.5508, -32.9514,2.77917, -0.0291596},{59.7514, -47.3033,3.54495, -0.0341802}}
            }
        };

        double[][][][] minparams = isinbending ? minparams_in : minparams_out;
        double[][][][] maxparams = isinbending ? maxparams_in : maxparams_out;


        double trajr = Math.sqrt(Math.pow(trajx,2) + Math.pow(trajy,2) + Math.pow(trajz,2));
        double theta_DCr = Math.toDegrees(Math.acos(trajz/trajr));
        double phi_DCr_raw = Math.toDegrees(Math.atan2(trajy/trajr, trajx/trajr));

        double phi_DCr = 5000;

        if (dc_sector == 1) phi_DCr = phi_DCr_raw;
        if (dc_sector == 2) phi_DCr = phi_DCr_raw - 60;
        if (dc_sector == 3) phi_DCr = phi_DCr_raw - 120;
        if (dc_sector == 4 && phi_DCr_raw > 0) phi_DCr = phi_DCr_raw - 180;
        if (dc_sector == 4 && phi_DCr_raw < 0) phi_DCr = phi_DCr_raw + 180;
        if (dc_sector == 5) phi_DCr = phi_DCr_raw + 120;
        if (dc_sector == 6) phi_DCr = phi_DCr_raw + 60;

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

        double calc_phi_min = minparams[pid][dc_sector-1][region-1][0] + minparams[pid][dc_sector-1][region-1][1] * Math.log(theta_DCr) + minparams[pid][dc_sector-1][region-1][2] * theta_DCr + minparams[pid][dc_sector-1][region-1][3] * theta_DCr * theta_DCr;

        double calc_phi_max = maxparams[pid][dc_sector-1][region-1][0] + maxparams[pid][dc_sector-1][region-1][1] * Math.log(theta_DCr) + maxparams[pid][dc_sector-1][region-1][2] * theta_DCr + maxparams[pid][dc_sector-1][region-1][3] * theta_DCr * theta_DCr;

        return (phi_DCr > calc_phi_min) && (phi_DCr < calc_phi_max);
    }



    /** Delta VZ cut for hadrons
     * @param pid hadron PID code
     * @param dvz difference between Vz of hadron candidate and electron
    */
    public static boolean Delta_vz_cut(int pid, double dvz) {
        switch(pid) {
        case 2212:
            return dvz>-20 && dvz<20;
        case 22:
            return dvz>-20 && dvz<20;
        case 2112:
            return dvz>-20 && dvz<20;
        case 211:
            return dvz>-20 && dvz<20;
        case -211:
            return dvz>-20 && dvz<20;
        case 321:
            return dvz>-20 && dvz<20;
        case -321:
            return dvz>-20 && dvz<20;
        }
        return false;
    }


}
