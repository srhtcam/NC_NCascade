using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;
using Rhino.Input.Custom;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Forms.VisualStyles;

namespace NC_NCascade
{
    public class NC4_NCascade : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public NC4_NCascade()
          : base("NC4_NCascade", "NC4_NCascade",
            "NC4_NCascade",
            "PSplines", "NCNCascade")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("pts", "pts", "pts", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("dij", "dij", "dij", GH_ParamAccess.list);
            pManager.AddLineParameter("dij_Lines", "dij_Lines", "dij_Lines", GH_ParamAccess.list);
            pManager.AddLineParameter("tensorLines", "tensorLines", "tensorLines", GH_ParamAccess.list);
            pManager.AddLineParameter("rptensorLines", "rptensorLines", "rptensorLines", GH_ParamAccess.list);
            pManager.AddPointParameter("repTBpts", "repTBpts", "repTBpts", GH_ParamAccess.list);
            pManager.AddGeometryParameter("nc_Pts", "nc_Pts", "nc_Pts", GH_ParamAccess.list);
            pManager.AddGeometryParameter("nc_Edges", "nc_Edges", "nc_Edges", GH_ParamAccess.list);
            pManager.AddGeometryParameter("nc_Surfs", "nc_Surfs", "nc_Surfs", GH_ParamAccess.list);
            pManager.AddTextParameter("bvfile", "bvfile", "bvfile", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Point3d> points = new List<Point3d>();
            if(!DA.GetDataList(0, points)) return;


            // Estimating n-Cascade size from input pts
            int a = points.Count, aa = 0;
            bool cas_size_find = true;
            int cas_size = 1;
            a = a - 14;
            while (a != 0)
            {
                cas_size++;
                a = a - (6 + aa);
                if (a == 0)break;
                else if (a < 0)
                {
                    cas_size_find = false;
                    break;
                }
                aa++;
            }

            if (!cas_size_find)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Net size does not match any n-Cascade size. Please check it. For referance;" +
                    "\n 1-cas net size:14," +
                    "\n 2-cas net size:20," +
                    "\n 3-cas net size:27," +
                    "\n 4-cas net size:35," +
                    "\n 5-cas net size:44 should be..");
                return;
            }

            // Assigning the input pts with respect to their cascade size
            int nn = cas_size + 3;
            Point3d[,] d_ij = new Point3d[nn, nn];
            for (int i = 0, ii = 0, mm = 0; i < nn; i++)
            {
                for (int j = 0; j < nn - mm; j++)
                {
                    d_ij[j, i] = points[ii];
                    ii++;
                }
                for (int j = nn - mm; j < nn; j++)
                {
                    d_ij[j, i] = Point3d.Unset;
                }
                if (i > 0 && i < nn - 2) mm++;
            }

            // Generating net from input pts.
            SetNet(d_ij,out List<Line> netlines); 

            // Constructing Tensor Border
            SetTensorBorder(d_ij, out List<Point3d[,]> t_ij, out List<Point3d[,]> tT_ij, out List<Point3d[,]> tB_ij, out List<Line> tblines);

            // Constructing Reparametrized Tensor Border
            SetReparametizedTensor(t_ij, tT_ij, tB_ij, out List<Point3d[,]> rp_t_ij, out List<Point3d[,]> rp_tB_ij, out List<Point3d[,]> rp_tT_ij, out List<Line> rp_tblines);

            // Reparametrized Tensor Border Points arrangenment for surface calcualtion, for using Matlab or another solver
            List<Point3d> rep_Tpts = RepTensorPtOut(rp_t_ij, rp_tB_ij, rp_tT_ij);

            // Arranged Rep. Tensor Border Points Reading
            Point3d[,] CP = ReptensorPtIn(rep_Tpts,cas_size);

            // Contiunity rules; BB-coefficients were calculated from rep. tensor border points
            Point3d[,] CP_final = ContinuityRulesNCasTC4(d_ij, CP);

            // Arranging the points for next surface construction method
            List<Point3d> outpts = new List<Point3d>();
            int m = 9 + (cas_size - 1) * 3;
            int n = (2 * (cas_size + 1)) + 1;
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    outpts.Add(CP_final[i, j]);
                }
            }


            // Surface construction from calculated points
            NC_nSurface(outpts, cas_size, out List<NurbsSurface> ncsurfs, out string bvfile, out List<Line> edges);

            DA.SetDataList(0, DAToList(d_ij));
            DA.SetDataList(1, netlines);
            DA.SetDataList(2, tblines);
            DA.SetDataList(3, rp_tblines);
            DA.SetDataList(4, rep_Tpts);
            DA.SetDataList(5, outpts);
            DA.SetDataList(6, edges);
            DA.SetDataList(7, ncsurfs);
            DA.SetData(8, bvfile);
        }

        private void NC_nSurface(List<Point3d> points,int n, out List<NurbsSurface> surfs, out string bvfile, out List<Line> tensoredges)
        {
            bvfile = "Group 2 NC4_"+n.ToString()+"Cas_Cap\n";
            surfs = new List<NurbsSurface>();

            int v_size = 3;
            int v_piece = n+1;
            int u_piece = n+1;

            int uu = 9 + ((u_piece - 2) * 3);
            int vv = ((v_size - 1) * v_piece) + 1;

            for (int ii = 0; ii < (vv - 2); ii = ii + 2)
            {
                int i = ii;
                surfs.Add(NurbsSurface.CreateFromPoints(new Point3d[] {
                    points[i], points[i + 1], points[i + 2],
                    points[i + vv], points[i + (vv + 1)], points[i + (vv + 2)],
                    points[i + 2*vv], points[i + (2*vv + 1)], points[i + (2*vv + 2)],
                    points[i + 3*vv], points[i + (3*vv + 1)], points[i + (3*vv + 2)],
                    points[i + 4*vv], points[i + (4*vv + 1)], points[i + (4*vv + 2)]}, 5, 3, 4, 2));

                bvfile += "5\n";
                bvfile += "4 2\n";
                bvfile += points[i].X + " " + points[i].Y + " " + points[i].Z + "\n";
                bvfile += points[i + 1].X + " " + points[i + 1].Y + " " + points[i + 1].Z + "\n";
                bvfile += points[i + 2].X + " " + points[i + 2].Y + " " + points[i + 2].Z + "\n";
                bvfile += points[i + vv].X + " " + points[i + vv].Y + " " + points[i + vv].Z + "\n";
                bvfile += points[i + (vv + 1)].X + " " + points[i + (vv + 1)].Y + " " + points[i + (vv + 1)].Z + "\n";
                bvfile += points[i + (vv + 2)].X + " " + points[i + (vv + 2)].Y + " " + points[i + (vv + 2)].Z + "\n";
                bvfile += points[i + 2 * vv].X + " " + points[i + 2 * vv].Y + " " + points[i + 2 * vv].Z + "\n";
                bvfile += points[i + (2 * vv + 1)].X + " " + points[i + (2 * vv + 1)].Y + " " + points[i + (2 * vv + 1)].Z + "\n";
                bvfile += points[i + (2 * vv + 2)].X + " " + points[i + (2 * vv + 2)].Y + " " + points[i + (2 * vv + 2)].Z + "\n";
                bvfile += points[i + 3 * vv].X + " " + points[i + 3 * vv].Y + " " + points[i + 3 * vv].Z + "\n";
                bvfile += points[i + (3 * vv + 1)].X + " " + points[i + (3 * vv + 1)].Y + " " + points[i + (3 * vv + 1)].Z + "\n";
                bvfile += points[i + (3 * vv + 2)].X + " " + points[i + (3 * vv + 2)].Y + " " + points[i + (3 * vv + 2)].Z + "\n";
                bvfile += points[i + 4 * vv].X + " " + points[i + 4 * vv].Y + " " + points[i + 4 * vv].Z + "\n";
                bvfile += points[i + (4 * vv + 1)].X + " " + points[i + (4 * vv + 1)].Y + " " + points[i + (4 * vv + 1)].Z + "\n";
                bvfile += points[i + (4 * vv + 2)].X + " " + points[i + (4 * vv + 2)].Y + " " + points[i + (4 * vv + 2)].Z + "\n";
            }

            for (int jj = 0; jj < (u_piece - 2); jj++)
            {
                for (int ii = ((4 + (jj * 3)) * vv); ii < (((4 + (jj * 3)) * vv) + (vv - 2)); ii = ii + 2)
                {
                    int i = ii;
                    surfs.Add(NurbsSurface.CreateFromPoints(new Point3d[] {
                    points[i], points[i + 1], points[i + 2],
                    points[i + vv], points[i + (vv + 1)], points[i + (vv + 2)],
                    points[i + 2*vv], points[i + (2*vv + 1)], points[i + (2*vv + 2)],
                    points[i + 3*vv], points[i + (3*vv + 1)], points[i + (3*vv + 2)] }, 4, 3, 3, 2));

                    bvfile += "5\n";
                    bvfile += "3 2\n";
                    bvfile += points[i].X + " " + points[i].Y + " " + points[i].Z + "\n";
                    bvfile += points[i + 1].X + " " + points[i + 1].Y + " " + points[i + 1].Z + "\n";
                    bvfile += points[i + 2].X + " " + points[i + 2].Y + " " + points[i + 2].Z + "\n";
                    bvfile += points[i + vv].X + " " + points[i + vv].Y + " " + points[i + vv].Z + "\n";
                    bvfile += points[i + (vv + 1)].X + " " + points[i + (vv + 1)].Y + " " + points[i + (vv + 1)].Z + "\n";
                    bvfile += points[i + (vv + 2)].X + " " + points[i + (vv + 2)].Y + " " + points[i + (vv + 2)].Z + "\n";
                    bvfile += points[i + 2 * vv].X + " " + points[i + 2 * vv].Y + " " + points[i + 2 * vv].Z + "\n";
                    bvfile += points[i + (2 * vv + 1)].X + " " + points[i + (2 * vv + 1)].Y + " " + points[i + (2 * vv + 1)].Z + "\n";
                    bvfile += points[i + (2 * vv + 2)].X + " " + points[i + (2 * vv + 2)].Y + " " + points[i + (2 * vv + 2)].Z + "\n";
                    bvfile += points[i + 3 * vv].X + " " + points[i + 3 * vv].Y + " " + points[i + 3 * vv].Z + "\n";
                    bvfile += points[i + (3 * vv + 1)].X + " " + points[i + (3 * vv + 1)].Y + " " + points[i + (3 * vv + 1)].Z + "\n";
                    bvfile += points[i + (3 * vv + 2)].X + " " + points[i + (3 * vv + 2)].Y + " " + points[i + (3 * vv + 2)].Z + "\n";
                }
            }


            for (int ii = ((uu - 5) * vv); ii < (((uu - 5) * vv) + (vv - 2)); ii = ii + 2)
            {
                int i = ii;
                surfs.Add(NurbsSurface.CreateFromPoints(new Point3d[] {
                    points[i], points[i + 1], points[i + 2],
                    points[i + vv], points[i + (vv + 1)], points[i + (vv + 2)],
                    points[i + 2*vv], points[i + (2*vv + 1)], points[i + (2*vv + 2)],
                    points[i + 3*vv], points[i + (3*vv + 1)], points[i + (3*vv + 2)],
                    points[i + 4*vv], points[i + (4*vv + 1)], points[i + (4*vv + 2)]}, 5, 3, 4, 2));

                bvfile += "5\n";
                bvfile += "4 2\n";
                bvfile += points[i].X + " " + points[i].Y + " " + points[i].Z + "\n";
                bvfile += points[i + 1].X + " " + points[i + 1].Y + " " + points[i + 1].Z + "\n";
                bvfile += points[i + 2].X + " " + points[i + 2].Y + " " + points[i + 2].Z + "\n";
                bvfile += points[i + vv].X + " " + points[i + vv].Y + " " + points[i + vv].Z + "\n";
                bvfile += points[i + (vv + 1)].X + " " + points[i + (vv + 1)].Y + " " + points[i + (vv + 1)].Z + "\n";
                bvfile += points[i + (vv + 2)].X + " " + points[i + (vv + 2)].Y + " " + points[i + (vv + 2)].Z + "\n";
                bvfile += points[i + 2 * vv].X + " " + points[i + 2 * vv].Y + " " + points[i + 2 * vv].Z + "\n";
                bvfile += points[i + (2 * vv + 1)].X + " " + points[i + (2 * vv + 1)].Y + " " + points[i + (2 * vv + 1)].Z + "\n";
                bvfile += points[i + (2 * vv + 2)].X + " " + points[i + (2 * vv + 2)].Y + " " + points[i + (2 * vv + 2)].Z + "\n";
                bvfile += points[i + 3 * vv].X + " " + points[i + 3 * vv].Y + " " + points[i + 3 * vv].Z + "\n";
                bvfile += points[i + (3 * vv + 1)].X + " " + points[i + (3 * vv + 1)].Y + " " + points[i + (3 * vv + 1)].Z + "\n";
                bvfile += points[i + (3 * vv + 2)].X + " " + points[i + (3 * vv + 2)].Y + " " + points[i + (3 * vv + 2)].Z + "\n";
                bvfile += points[i + 4 * vv].X + " " + points[i + 4 * vv].Y + " " + points[i + 4 * vv].Z + "\n";
                bvfile += points[i + (4 * vv + 1)].X + " " + points[i + (4 * vv + 1)].Y + " " + points[i + (4 * vv + 1)].Z + "\n";
                bvfile += points[i + (4 * vv + 2)].X + " " + points[i + (4 * vv + 2)].Y + " " + points[i + (4 * vv + 2)].Z + "\n";
            }

            tensoredges = new List<Line>();
            // Horizontal
            for (int i = 0; i < (uu * vv); i++)
            {
                if (i % vv == (vv - 1)) continue;
                tensoredges.Add(new Line(points[i], points[i + 1]));
            }

            // Vertical
            for (int i = 0; i < ((uu - 1) * vv); i++)
            {
                tensoredges.Add(new Line(points[i], points[i + vv]));
            }
        }

        public static Point3d[,] ReptensorPtIn(List<Point3d> pts,int cas_n)
        {
            cas_n++;// npc size ..

            int m = (9)+((cas_n - 2)*3), n = 2*cas_n +1;
            Point3d[,] CP = new Point3d[m, n];

            int indx = 0;
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < m; i++)
                {
                    if (j < 2)
                    {
                        CP[i, j] = pts[indx];
                        indx++;
                    }
                    else if (j > 1 && j < (n - 2))
                    {
                        CP[0, j] = pts[indx]; indx++;
                        CP[1, j] = pts[indx]; indx++;
                        CP[(m - 2), j] = pts[indx]; indx++;
                        CP[(m - 1), j] = pts[indx]; indx++;
                        break;
                    }
                    else if (j > (n - 3))
                    {
                        CP[i, j] = pts[indx];
                        indx++;
                    }
                }

            }
            return CP;
        }

        public static Point3d[,] ContinuityRulesNCasTC4(Point3d[,] netpts, Point3d[,] CP_in)
        {
            int u = CP_in.GetLength(0);
            int v = CP_in.GetLength(1);
            Point3d[,] CP = new Point3d[u, v];

            for (int i = 0; i < u; i++)
            {
                for (int j = 0; j < v; j++)
                {
                    CP[i, j] = new Point3d(CP_in[i, j]); // Copy each input element to new list
                }
            }
            
            ///
            int n = (int)Math.Floor(v/2.0);

            List<List<Vector3d>> c = new List<List<Vector3d>>();

            for (int i = 0; i < n - 2; i++)
            {
                c.Add(new List<Vector3d>());
                for (int j = 0; j < n; j++)
                {
                    if (j == 0)
                        c[c.Count - 1].Add((2 * CP[1, (2 * i + 3)]) - CP[0, (2 * i + 3)]);
                    else if (j == n - 1)
                        c[c.Count - 1].Add((2 * CP[(u - 2), (2 * i + 3)]) - CP[(u - 1), (2 * i + 3)]);
                    else
                        c[c.Count - 1].Add(CalculateCij_kk(i+2, j+1, n, netpts));
                }
            }

            int jj = 0;
            for (int j = 3; j < v - 2; j = j + 2)
            {
                // Diamond rule and bottom bi-2
                CP[4, j] = (Point3d)(0.5 * (c[jj][0] + c[jj][1]));
                CP[3, j] = 0.5 * (c[jj][0] + CP[4, j]);
                CP[2, j] = (1.0 / 6) * CP[0, j] + (2.0 / 3) * c[jj][0] + (1.0 / 6) * CP[4, j];

                // Diamond rule and top bi-2
                CP[(u - 5), j] = (Point3d)(0.5 * (c[jj][n - 1] + c[jj][n - 2]));
                CP[(u - 4), j] = 0.5 * (c[jj][n - 1] + CP[(u - 5), j]);
                CP[(u - 3), j] = (1.0 / 6) * CP[(u - 1), j] + (2.0 / 3) * c[jj][n - 1] + (1.0 / 6) * CP[(u - 5), j];

                // middle patch vertical spine pieces
                int ii = 1;
                for (int i = 7; i <= (u - 5); i = i + 3)
                {
                    CP[i, j] = (Point3d)(0.5 * (c[jj][ii] + c[jj][ii + 1]));

                    CP[(i - 2), j] = (1.0 / 3) * CP[(i - 3), j] + (2.0 / 3) * c[jj][ii];
                    CP[(i - 1), j] = (2.0 / 3) * c[jj][ii] + (1.0 / 3) * CP[i, j];

                    ii++;// vertical spine piece increament for "c" list
                }
                jj++;// next spine increament for "c" list
            }


            // g1 contiunity
            for (int i = 2; i <= u - 3; i++)
            {
                for (int j = 2; j < v - 2; j = j + 2)
                {
                    CP[i, j] = 0.5 * (CP[i, (j - 1)] + CP[i, (j + 1)]);
                }
            }

            return CP;
        }

        private static Vector3d CalculateCij_kk(int i, int j, int n, Point3d[,] net)
        {
            double r = (((n + 1 - j) / (2.0 * n)) * (2 * i - 1)) + 1;
            int s = (int)Math.Floor(r);
            double u = r - s;

            Vector3d cij = (Vector3d)(((net[s - 1, j] + net[s, j]) / 2) * (Math.Pow((1 - u), 2)) +
                   (net[s, j]) * (2 * u * (1 - u)) +
                   ((net[s, j] + net[s + 1, j]) / 2) * (Math.Pow(u, 2)));

            return cij;
        }

        public void SetReparametizedTensor(List<Point3d[,]> t_ij, List<Point3d[,]> tT_ij, List<Point3d[,]> tB_ij, out List<Point3d[,]> rp_t_ij, out List<Point3d[,]> rp_tB_ij, out List<Point3d[,]> rp_tT_ij, out List<Line> rp_tblines)
        {
            rp_t_ij = new List<Point3d[,]>();
            rp_tB_ij = new List<Point3d[,]>();
            rp_tT_ij = new List<Point3d[,]>();

            #region left and right

            int n = t_ij.Count/2;

            List<double> a0 = new List<double>();
            List<double> a1 = new List<double>();
            List<double> a2 = new List<double>();

            a0.Add(1.0);
            a1.Add(1.0);
            a2.Add(1.0 - ((1.0) / (2 * n)));
            for (int i = 2; i < n; i++)
            {
                a0.Add(1.0 - (((2 * i) - 3.0) / (2 * n)));
                a1.Add(1.0 - (((2 * (i + 1)) - 3.0) / (2 * n)));
                a2.Add(0);
            }
            a0.Add(1.0 - (((2 * n) - 3.0) / (2 * n)));
            a1.Add(1.0 / n);
            a2.Add(1.0 / n);

            //right tensor border coefficients
            a0.AddRange(a0);
            a1.AddRange(a1);
            a2.AddRange(a2);

            for (int i = 0; i < 2*n; i++)
            {
                Point3d par_t00, par_t10, par_t20, par_t30, par_t40, par_t01, par_t11, par_t21, par_t31, par_t41;
                Point3d
                    t00 = t_ij[i][0, 0],
                    t10 = t_ij[i][0, 1],
                    t20 = t_ij[i][0, 2],
                    t01 = t_ij[i][1, 0],
                    t11 = t_ij[i][1, 1],
                    t21 = t_ij[i][1, 2];

                if (i > n-1) // right tensor border
                {
                    t00 = t_ij[i][1, 0];
                    t10 = t_ij[i][1, 1];
                    t20 = t_ij[i][1, 2];
                    t01 = t_ij[i][0, 0];
                    t11 = t_ij[i][0, 1];
                    t21 = t_ij[i][0, 2];
                }

                if (i == 0 || i == (n - 1) || i == n || i == (2 * n - 1))
                {
                    par_t00 = t00;
                    par_t10 = t00 / 2 + t10 / 2;
                    par_t20 = t00 / 6 + 2 * t10 / 3 + t20 / 6;
                    par_t30 = t10 / 2 + t20 / 2;
                    par_t40 = t20;
                    par_t01 = -a0[i] * t00 + a0[i] * t01 + t00;
                    par_t11 = -a0[i] * t10 / 2 + a0[i] * t11 / 2 - a1[i] * t00 / 2 + a1[i] * t01 / 2 + t00 / 2 + t10 / 2;
                    par_t21 = -a0[i] * t20 / 6 + a0[i] * t21 / 6 - 2 * a1[i] * t10 / 3 + 2 * a1[i] * t11 / 3 - a2[i] * t00 / 6 + a2[i] * t01 / 6 + t00 / 6 + 2 * t10 / 3 + t20 / 6;
                    par_t31 = -a1[i] * t20 / 2 + a1[i] * t21 / 2 - a2[i] * t10 / 2 + a2[i] * t11 / 2 + t10 / 2 + t20 / 2;
                    par_t41 = -a2[i] * t20 + a2[i] * t21 + t20;

                    rp_t_ij.Add(new Point3d[2, 5]);
                    rp_t_ij[i][0, 0] = par_t00;
                    rp_t_ij[i][0, 1] = par_t10;
                    rp_t_ij[i][0, 2] = par_t20;
                    rp_t_ij[i][0, 3] = par_t30;
                    rp_t_ij[i][0, 4] = par_t40;
                    rp_t_ij[i][1, 0] = par_t01;
                    rp_t_ij[i][1, 1] = par_t11;
                    rp_t_ij[i][1, 2] = par_t21;
                    rp_t_ij[i][1, 3] = par_t31;
                    rp_t_ij[i][1, 4] = par_t41;
                }
                else
                {
                    par_t00 = t00;
                    par_t10 = t00 / 3 + 2 * t10 / 3;
                    par_t20 = 2 * t10 / 3 + t20 / 3;
                    par_t30 = t20;
                    par_t01 = -a0[i] * t00 + a0[i] * t01 + t00;
                    par_t11 = -2 * a0[i] * t10 / 3 + 2 * a0[i] * t11 / 3 - a1[i] * t00 / 3 + a1[i] * t01 / 3 + t00 / 3 + 2 * t10 / 3;
                    par_t21 = -a0[i] * t20 / 3 + a0[i] * t21 / 3 - 2 * a1[i] * t10 / 3 + 2 * a1[i] * t11 / 3 + 2 * t10 / 3 + t20 / 3;
                    par_t31 = -a1[i] * t20 + a1[i] * t21 + t20;

                    rp_t_ij.Add(new Point3d[2, 4]);
                    rp_t_ij[i][0, 0] = par_t00;
                    rp_t_ij[i][0, 1] = par_t10;
                    rp_t_ij[i][0, 2] = par_t20;
                    rp_t_ij[i][0, 3] = par_t30;
                    rp_t_ij[i][1, 0] = par_t01;
                    rp_t_ij[i][1, 1] = par_t11;
                    rp_t_ij[i][1, 2] = par_t21;
                    rp_t_ij[i][1, 3] = par_t31;
                }

            }

            //rep. tensor line
            rp_tblines = new List<Line>();
            for (int i = 0; i < rp_t_ij.Count; i++)
            {
                for (int j = 0; j < rp_t_ij[i].GetLength(1); j++)
                {
                    rp_tblines.Add(new Line(rp_t_ij[i][0, j], rp_t_ij[i][1, j]));
                }
                rp_tblines.Add(new Line(rp_t_ij[i][0, 0], rp_t_ij[i][0, 1]));
                rp_tblines.Add(new Line(rp_t_ij[i][0, 1], rp_t_ij[i][0, 2]));
                rp_tblines.Add(new Line(rp_t_ij[i][1, 0], rp_t_ij[i][1, 1]));
                rp_tblines.Add(new Line(rp_t_ij[i][1, 1], rp_t_ij[i][1, 2]));
                rp_tblines.Add(new Line(rp_t_ij[i][0, 2], rp_t_ij[i][0, 3]));
                rp_tblines.Add(new Line(rp_t_ij[i][1, 2], rp_t_ij[i][1, 3]));
                if(i == 0 || i == (n - 1) || i == n || i == (2 * n - 1))
                {
                    rp_tblines.Add(new Line(rp_t_ij[i][0, 3], rp_t_ij[i][0, 4]));
                    rp_tblines.Add(new Line(rp_t_ij[i][1, 3], rp_t_ij[i][1, 4]));
                }
            }

            #endregion

            #region tB
            ///Bt
            ///
            /////
            for (int i = 0; i < n; i++)
            {
                rp_tB_ij.Add(new Point3d[3, 2]);
                rp_tB_ij[i][0, 0] = tB_ij[i][0, 0];
                rp_tB_ij[i][1, 0] = tB_ij[i][1, 0];
                rp_tB_ij[i][2, 0] = tB_ij[i][2, 0];
                rp_tB_ij[i][0, 1] = (tB_ij[i][0, 0] / 2) + (tB_ij[i][0, 1] / 2);
                rp_tB_ij[i][1, 1] = (tB_ij[i][1, 0] / 2) + (tB_ij[i][1, 1] / 2);
                rp_tB_ij[i][2, 1] = (tB_ij[i][2, 0] / 2) + (tB_ij[i][2, 1] / 2);
            }

            //rep. tensor line
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    rp_tblines.Add(new Line(rp_tB_ij[i][j, 0], rp_tB_ij[i][j, 1]));
                }
                for (int j = 0; j < 2; j++)
                {
                    rp_tblines.Add(new Line(rp_tB_ij[i][0, j], rp_tB_ij[i][1, j]));
                    rp_tblines.Add(new Line(rp_tB_ij[i][1, j], rp_tB_ij[i][2, j]));
                }
            }
            #endregion

            #region tT
            ///tT
            /////

            //bottom
            List<Point3d[]> newcntrolpts_bottom = new List<Point3d[]>();
            List<Point3d[]> decastel_bottom = new List<Point3d[]>();
            decastel_bottom.Add(new Point3d[] { tT_ij[0][0, 0], tT_ij[0][1, 0], tT_ij[0][2, 0] });

            //top
            List<Point3d[]> newcntrolpts_top = new List<Point3d[]>();
            List<Point3d[]> decastel_top = new List<Point3d[]>();
            decastel_top.Add(new Point3d[] { tT_ij[0][0, 1], tT_ij[0][1, 1], tT_ij[0][2, 1] });

            //decasteljau
            for (int i = n; i > 1; i--)
            {
                DecasteljauDivision(decastel_bottom[decastel_bottom.Count - 1], (1.0 / i), out Point3d[] newcntrlpt0_bottom, out Point3d[] newcntrlpt1_bottom);
                decastel_bottom.Add(newcntrlpt1_bottom);
                newcntrolpts_bottom.Add(newcntrlpt0_bottom);

                DecasteljauDivision(decastel_top[decastel_top.Count - 1], (1.0 / i), out Point3d[] newcntrlpt0_top, out Point3d[] newcntrlpt1_top);
                decastel_top.Add(newcntrlpt1_top);
                newcntrolpts_top.Add(newcntrlpt0_top);
            }
            newcntrolpts_bottom.Add(decastel_bottom[decastel_bottom.Count - 1]);
            newcntrolpts_top.Add(decastel_top[decastel_top.Count - 1]);


            for (int i = 0; i < n; i++)
            {
                rp_tT_ij.Add(new Point3d[3, 2]);
                rp_tT_ij[i][0, 0] = (newcntrolpts_bottom[i][0] / 2) + (newcntrolpts_top[i][0] / 2);
                rp_tT_ij[i][1, 0] = (newcntrolpts_bottom[i][1] / 2) + (newcntrolpts_top[i][1] / 2);
                rp_tT_ij[i][2, 0] = (newcntrolpts_bottom[i][2] / 2) + (newcntrolpts_top[i][2] / 2);
                rp_tT_ij[i][0, 1] = newcntrolpts_top[i][0];
                rp_tT_ij[i][1, 1] = newcntrolpts_top[i][1];
                rp_tT_ij[i][2, 1] = newcntrolpts_top[i][2];
            }

            //rep. tensor line
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    rp_tblines.Add(new Line(rp_tT_ij[i][j, 0], rp_tT_ij[i][j, 1]));
                }
                for (int j = 0; j < 2; j++)
                {
                    rp_tblines.Add(new Line(rp_tT_ij[i][0, j], rp_tT_ij[i][1, j]));
                    rp_tblines.Add(new Line(rp_tT_ij[i][1, j], rp_tT_ij[i][2, j]));
                }
            }
            ///
            #endregion
        }

        private List<Point3d> RepTensorPtOut(List<Point3d[,]> rp_t_ij, List<Point3d[,]> rp_tB_ij, List<Point3d[,]> rp_tT_ij)
        {
            List<Point3d> pts = new List<Point3d>();
            int n = rp_tB_ij.Count;

            for (int i = 0; i < n; i++)// vertical line 0
            {
                if (i == 0) pts.Add(rp_t_ij[i][0, 0]);
                pts.Add(rp_t_ij[i][0, 1]);
                pts.Add(rp_t_ij[i][0, 2]);
                pts.Add(rp_t_ij[i][0, 3]);
                if (i == 0 || i == n-1) pts.Add(rp_t_ij[i][0, 4]);
            }
            for (int i = 0; i < n; i++)//vertical line 1
            {
                if (i == 0) pts.Add(rp_t_ij[i][1, 0]);
                pts.Add(rp_t_ij[i][1, 1]);
                pts.Add(rp_t_ij[i][1, 2]);
                pts.Add(rp_t_ij[i][1, 3]);
                if (i == 0 || i == n-1) pts.Add(rp_t_ij[i][1, 4]);
            }
            for (int i = 0; i < n-1; i++)
            {
                if (i != 0)
                {
                    pts.Add(rp_tB_ij[i][1, 0]);
                    pts.Add(rp_tB_ij[i][1, 1]);
                    pts.Add(rp_tT_ij[i][1, 0]);
                    pts.Add(rp_tT_ij[i][1, 1]);
                }
                if (i != n-1)
                {
                    pts.Add(rp_tB_ij[i][2, 0]);
                    pts.Add(rp_tB_ij[i][2, 1]);
                    pts.Add(rp_tT_ij[i][2, 0]);
                    pts.Add(rp_tT_ij[i][2, 1]);
                }

            }
            for (int i = n; i < 2*n; i++)//vertical line 2n-1
            {
                if (i == n) pts.Add(rp_t_ij[i][1, 0]);
                pts.Add(rp_t_ij[i][1, 1]);
                pts.Add(rp_t_ij[i][1, 2]);
                pts.Add(rp_t_ij[i][1, 3]);
                if (i == n || i == 2*n-1) pts.Add(rp_t_ij[i][1, 4]);
            }
            for (int i = n; i < 2*n; i++) //vertical line 2n
            {
                if (i == n) pts.Add(rp_t_ij[i][0, 0]);
                pts.Add(rp_t_ij[i][0, 1]);
                pts.Add(rp_t_ij[i][0, 2]);
                pts.Add(rp_t_ij[i][0, 3]);
                if (i == n || i == 2*n-1) pts.Add(rp_t_ij[i][0, 4]);
            }
            return pts;
        }

        /// <summary>
        /// https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/de-casteljau.html
        /// </summary>
        /// <param name="controlpt"></param>
        /// <param name="u"></param>
        /// <param name="newcontrolptUpper"></param>
        /// <param name="newcontrolptBottom"></param>
        public void DecasteljauDivision(Point3d[] controlpt, double u, out Point3d[] newcontrolptUpper, out Point3d[] newcontrolptBottom)
        {
            newcontrolptUpper = new Point3d[controlpt.Length];
            newcontrolptBottom = new Point3d[controlpt.Length];

            List<List<Point3d>> ptMatrix = new List<List<Point3d>>();
            ptMatrix.Add(new List<Point3d>(controlpt.ToList()));
            for (int i = 0; i < controlpt.Length; i++)
            {
                ptMatrix.Add(new List<Point3d>());
                for (int j = 0; j < ptMatrix[i].Count - 1; j++)
                {
                    ptMatrix[i + 1].Add((1 - u) * ptMatrix[i][j] + u * ptMatrix[i][j + 1]);
                }
                newcontrolptUpper[i] = new Point3d(ptMatrix[i][0]);
                newcontrolptBottom[(controlpt.Length - 1) - i] = new Point3d(ptMatrix[i][ptMatrix[i].Count - 1]);
            }
        }

        private void SetTensorBorder(Point3d[,] d_ij, out List<Point3d[,]> t_ij, out List<Point3d[,]> tT_ij, out List<Point3d[,]> tB_ij, out List<Line> tblines)
        {
            t_ij = new List<Point3d[,]>();
            tT_ij = new List<Point3d[,]>();
            tB_ij = new List<Point3d[,]>();

            int n = d_ij.GetLength(0) - 2;

            //left
            for (int i = 0; i < n; i++)
            {
                t_ij.Add(new Point3d[2, 3]);
                t_ij[i][0, 0] = SetQMiddle(d_ij[0, i + 0], d_ij[0, i + 1], d_ij[1, i + 0], d_ij[1, i + 1]);
                t_ij[i][1, 0] = SetEMiddle(d_ij[1, i + 0], d_ij[1, i + 1]);
                t_ij[i][0, 1] = SetEMiddle(d_ij[0, i + 1], d_ij[1, i + 1]);
                t_ij[i][1, 1] = d_ij[1, i + 1];
                t_ij[i][0, 2] = SetQMiddle(d_ij[0, i + 1], d_ij[0, i + 2], d_ij[1, i + 1], d_ij[1, i + 2]);
                t_ij[i][1, 2] = SetEMiddle(d_ij[1, i + 1], d_ij[1, i + 2]);
            }

            //right
            t_ij.Add(new Point3d[2, 3]);
            t_ij[t_ij.Count - 1][0, 0] = SetEMiddle(d_ij[n, 0], d_ij[n, 1]);
            t_ij[t_ij.Count - 1][1, 0] = SetQMiddle(d_ij[n, 0], d_ij[n, 1], d_ij[n+1, 0], d_ij[n+1, 1]);
            t_ij[t_ij.Count - 1][0, 1] = d_ij[n, 1];
            t_ij[t_ij.Count - 1][1, 1] = SetEMiddle(d_ij[n, 1], d_ij[n+1, 1]);
            t_ij[t_ij.Count - 1][0, 2] = SetEMiddle(d_ij[n, 1], d_ij[n-1, 2]);
            t_ij[t_ij.Count - 1][1, 2] = SetQMiddle(d_ij[n, 1], d_ij[n-1, 2], d_ij[n+1, 1], d_ij[n, 2]);

            for (int i = 0; i < n - 2; i++)
            {
                t_ij.Add(new Point3d[2, 3]);
                t_ij[t_ij.Count - 1][0, 0] = SetEMiddle(d_ij[(n - i), 1 + i], d_ij[(n - 1) - i, 2 + i]);
                t_ij[t_ij.Count - 1][1, 0] = SetQMiddle(d_ij[(n - i), 1 + i], d_ij[(n - 1) - i, 2 + i], d_ij[(n + 1) - i, 1 + i], d_ij[(n) - i, 2 + i]);
                t_ij[t_ij.Count - 1][0, 1] = d_ij[(n - 1) - i, 2 + i];
                t_ij[t_ij.Count - 1][1, 1] = SetEMiddle(d_ij[(n - 1) - i, 2 + i], d_ij[(n) - i, 2 + i]);
                t_ij[t_ij.Count - 1][0, 2] = SetEMiddle(d_ij[(n - 1) - i, 2 + i], d_ij[(n - 2) - i, 3 + i]);
                t_ij[t_ij.Count - 1][1, 2] = SetQMiddle(d_ij[(n - 1) - i, 2 + i], d_ij[(n - 2) - i, 3 + i], d_ij[n - i, 2 + i], d_ij[(n - 1) - i, 3 + i]);
            }

            t_ij.Add(new Point3d[2, 3]);
            t_ij[t_ij.Count - 1][0, 0] = SetEMiddle(d_ij[2, n - 1], d_ij[1, n]);
            t_ij[t_ij.Count - 1][1, 0] = SetQMiddle(d_ij[2, n - 1], d_ij[1, n], d_ij[3, n - 1], d_ij[2, n]);
            t_ij[t_ij.Count - 1][0, 1] = d_ij[1, n];
            t_ij[t_ij.Count - 1][1, 1] = SetEMiddle(d_ij[1, n], d_ij[2, n]);
            t_ij[t_ij.Count - 1][0, 2] = SetEMiddle(d_ij[1, n], d_ij[1, n + 1]);
            t_ij[t_ij.Count - 1][1, 2] = SetQMiddle(d_ij[1, n], d_ij[1, n + 1], d_ij[2, n], d_ij[2, n + 1]);

            //top
            tT_ij.Add(new Point3d[3, 3]);
            tT_ij[0][0, 0] = SetEMiddle(d_ij[0, n], d_ij[1, n]);
            tT_ij[0][1, 0] = d_ij[1, n];
            tT_ij[0][2, 0] = SetEMiddle(d_ij[1, n], d_ij[2, n]);
            tT_ij[0][0, 1] = SetQMiddle(d_ij[0, n], d_ij[0, n + 1], d_ij[1, n], d_ij[1, n + 1]);
            tT_ij[0][1, 1] = SetEMiddle(d_ij[1, n], d_ij[1, n + 1]);
            tT_ij[0][2, 1] = SetQMiddle(d_ij[1, n], d_ij[1, n + 1], d_ij[2, n], d_ij[2, n + 1]);

            //bottom
            for (int i = 0; i < n; i++)
            {
                tB_ij.Add(new Point3d[3, 2]);
                tB_ij[i][0, 0] = SetQMiddle(d_ij[0 + i, 0], d_ij[0 + i, 1], d_ij[1 + i, 0], d_ij[1 + i, 1]);
                tB_ij[i][1, 0] = SetEMiddle(d_ij[1 + i, 0], d_ij[1 + i, 1]);
                tB_ij[i][2, 0] = SetQMiddle(d_ij[1 + i, 0], d_ij[1 + i, 1], d_ij[2 + i, 0], d_ij[2 + i, 1]);
                tB_ij[i][0, 1] = SetEMiddle(d_ij[0 + i, 1], d_ij[1 + i, 1]);
                tB_ij[i][1, 1] = d_ij[1 + i, 1];
                tB_ij[i][2, 1] = SetEMiddle(d_ij[1 + i, 1], d_ij[2 + i, 1]);
            }

            //tblines
            tblines = new List<Line>();
            //left-right
            for (int i = 0; i < 2*n; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    tblines.Add(new Line(t_ij[i][0, j], t_ij[i][1, j]));
                }
                tblines.Add(new Line(t_ij[i][0, 0], t_ij[i][0, 1]));
                tblines.Add(new Line(t_ij[i][0, 1], t_ij[i][0, 2]));
                tblines.Add(new Line(t_ij[i][1, 0], t_ij[i][1, 1]));
                tblines.Add(new Line(t_ij[i][1, 1], t_ij[i][1, 2]));
            }
            //bottom
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    tblines.Add(new Line(tB_ij[i][j, 0], tB_ij[i][j, 1]));
                }
                for (int j = 0; j < 2; j++)
                {
                    tblines.Add(new Line(tB_ij[i][0, j], tB_ij[i][1, j]));
                    tblines.Add(new Line(tB_ij[i][1, j], tB_ij[i][2, j]));
                }
            }
            //top
            for (int j = 0; j < 3; j++)
            {
                tblines.Add(new Line(tT_ij[0][j, 0], tT_ij[0][j, 1]));
            }
            for (int j = 0; j < 2; j++)
            {
                tblines.Add(new Line(tT_ij[0][0, j], tT_ij[0][1, j]));
                tblines.Add(new Line(tT_ij[0][1, j], tT_ij[0][2, j]));
            }
        }

        private Point3d SetQMiddle(Point3d p0, Point3d p1, Point3d p2, Point3d p3)
        {
            return (p0 + p1 + p2 + p3) / 4;
        }
        private Point3d SetEMiddle(Point3d p0, Point3d p1)
        {
            return (p0 + p1) / 2;
        }

        private void SetNet(Point3d[,] d_ij,out List<Line> netline)
        {
            netline = new List<Line>();
            
            int m = d_ij.GetLength(0), mm = 1, mm1 = 0;
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < m - mm1; j++)// Vertical lines
                {
                    if (i == m - 1) break;
                    else if (i == 0) netline.Add(new Line(d_ij[j, i], d_ij[j, i + 1]));
                    else if (i == m - 2) netline.Add(new Line(d_ij[j, i], d_ij[j, i + 1]));
                    else
                    {
                        if (j < 2) netline.Add(new Line(d_ij[j, i], d_ij[j, i + 1]));
                        else if (j > (m - mm1) - 3) netline.Add(new Line(d_ij[j, i], d_ij[j - 1, i + 1]));
                        else
                        {
                            netline.Add(new Line(d_ij[j, i], d_ij[j - 1, i + 1]));
                            netline.Add(new Line(d_ij[j, i], d_ij[j, i + 1]));
                        }
                    }
                }
                for (int j = 0; j < m - mm; j++)// Horizontal lines
                {
                    netline.Add(new Line(d_ij[j, i], d_ij[j + 1, i]));
                }
                if (i > 1 && i < m - 2) mm++;
                if(i > 0) mm1++;
            }
        }

        public static List<Point3d> DAToList(Point3d[,] array)
        {
            int width = array.GetLength(0);
            int height = array.GetLength(1);
            List<Point3d> ret = new List<Point3d>(width * height);
            for (int i = 0; i < width; i++)
            {
                for (int j = 0; j < height; j++)
                {
                    ret.Add(array[i, j]);
                }
            }
            return ret;
        }


        public override void AddRuntimeMessage(GH_RuntimeMessageLevel level, string text)
        {
            base.AddRuntimeMessage(level, text);
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// You can add image files to your project resources and access them like this:
        /// return Resources.IconForThisComponent;
        /// </summary>
        protected override System.Drawing.Bitmap Icon => null;

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid => new Guid("abf811e0-c174-4761-8bb9-cc85e9905de9");
    }
}