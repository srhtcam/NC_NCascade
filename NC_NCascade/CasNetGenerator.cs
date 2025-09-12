using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using Rhino.Input.Custom;
using System.Text;

namespace NC_NCascade
{
    public class CasNetGenerator : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public CasNetGenerator()
          : base("CasNetGenerator", "CasNetGenerator",
              "CasNetGenerator",
              "PSplines", "NCNCascade")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddIntegerParameter("n_size", "n_size", "n_size", GH_ParamAccess.item, 3);
            pManager.AddIntegerParameter("case", "case", "case", GH_ParamAccess.item, 0);
            pManager.AddNumberParameter("sph_r_inc", "sph_r_inc", "Sphere radius increment or decrement", GH_ParamAccess.item, 0);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("dij", "dij", "dij", GH_ParamAccess.list);
            pManager.AddLineParameter("dij_edge", "dij_edge", "dij_edge", GH_ParamAccess.list);
            pManager.AddGeometryParameter("sphere", "sphere", "sphere", GH_ParamAccess.item);
            pManager.AddGeometryParameter("extSurf", "extSurf", "extSurf", GH_ParamAccess.list);
            pManager.AddTextParameter("bv_file", "bv_file", "bv_file", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            int caseno = 0, cas_size = 3;
            double sphere_inc = 0;
            if (!DA.GetData(0, ref cas_size) || cas_size < 1) return;
            if (!DA.GetData(1, ref caseno)) return;
            if (!DA.GetData(2, ref sphere_inc)) return;

            // Generating a cascade net for a given size
            SetCasNet(cas_size, caseno, sphere_inc, out List<Point3d> netpts, out List<Point3d> net_extent, out Sphere sphere);

            // Constructing cascade net edges
            SetExtNet(cas_size, net_extent, out List<Line> extNetLines);

            // Constructing extended cascade net surfaces
            SetExtendedNCasSurf(cas_size, net_extent, out List<Point3d> tensorpoints, out List<NurbsSurface> surfs, out StringBuilder bvfile);

            DA.SetDataList(0, netpts);
            DA.SetDataList(1, extNetLines);
            DA.SetData(2, sphere);
            DA.SetDataList(3, surfs);
            DA.SetData(4, bvfile.ToString());
        }

        private void SetCasNet(int cas_size, int caseNo,double sphere_inc, out List<Point3d> netpts, out List<Point3d> net_extent, out Sphere sph)
        {
            // Net generator....
            int n = cas_size + 3;
            net_extent = new List<Point3d>();
            netpts = new List<Point3d>();

            double jj = 0;
            for (int i = 0; i < (n+2); i++)
            {
                for (double j = 0+jj; j < (n+2)-jj; j++)
                {
                    net_extent.Add(new Point3d(j, i, 0));
                }
                if(i > 1 && i < (n + 2)-3) jj = jj + 0.5;
            }

            // Projecting the net points to the sphere
            double sphere_r = 2* (n + 1) + sphere_inc;
            sph = new Sphere(new Point3d((n + 1) / 2.0, (n + 1) / 2.0, -(sphere_r + 1)), sphere_r);
            Line ln = new Line();
            Point3d Pt = new Point3d(), intPt1, intPt2;
            for (int i = 0; i < net_extent.Count; i++)
            {
                ln.From = net_extent[i];
                Pt = net_extent[i];
                Pt.Z = -sphere_r / 2.0;
                ln.To = Pt;

                Rhino.Geometry.Intersect.Intersection.LineSphere(ln, sph, out intPt1, out intPt2);
                net_extent[i] = intPt1;
            }

            // deformation cases
            if (caseNo == 0)
            {
                netpts = NetPtsOut(net_extent, n);
            }
            else if (cas_size == 3)
            {
                if (caseNo == 1)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 0.5;

                    int[] idx = new int[] { 25, 26, 27, 28, 29, 24, 30 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);

                }
                else if (caseNo == 2)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 0.5;

                    int[] idx = new int[] { 32, 33, 27, 12, 20, 31, 4 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 3)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = -1;

                    int[] idx = new int[] { 9, 17, 25, 32, 38, 43, 1, 48, 0, 8, 16, 24, 31, 37, 42, 47 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 4)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = -1;

                    int[] idx = new int[] { 9, 17, 25, 32, 38, 43, 40, 45, 35, 29, 14, 22, 1, 48, 0, 8, 16, 24, 31, 37, 42, 47, 6, 50, 7, 15, 23, 30, 36, 41, 46, 51 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 5)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 0.5;

                    int[] idx = new int[] { 11, 19, 27, 34, 12, 20, 28, 13, 21, 3, 4, 5 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 6)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 0.5;

                    int[] idx = new int[] { 11, 19, 27, 12, 20, 3, 4 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 7)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 4;

                    int[] idx = new int[] { 27 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);

                }
                else if (caseNo == 71)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = -0.5;

                    int[] idx = new int[] { 27 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);

                }
                else if (caseNo == 8)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = -0.5;
                    double translatedist2 = 0.5;

                    int[] idx = new int[] { 17, 32, 18, 33, 19, 34, 20, 35, 21, 22, 16, 23, 31, 36 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 9)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = -0.5;

                    int[] idx = new int[] { 17, 32, 43, 18, 33, 44, 19, 34, 45, 20, 35, 21, 22, 16, 23, 31, 36, 42, 46 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 10)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = -1;

                    int[] idx = new int[] { 11, 19, 27, 12, 20, 28, 13, 21, 29, 14, 22, 3, 4, 5, 6, 7, 15, 23, 30 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 101)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = -0.5;

                    int[] idx = new int[] { 11, 19, 27, 12, 20, 28, 13, 21, 29, 14, 22, 3, 4, 5, 6, 7, 15, 23, 30 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 11)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 1;

                    int[] idx = new int[] { 33, 39, 44, 23, 34, 49 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 111)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = -0.5;

                    int[] idx = new int[] { 33, 39, 44, 27, 34, 49 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 12)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = -1;


                    int[] idx = new int[] { 25, 32, 26, 33, 11, 19, 27, 34, 12, 20, 28, 35, 29, 3, 4, 24, 31, 36, 30 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 13)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 1;

                    int[] idx = new int[] { 25, 26, 33, 39, 44, 27, 34, 28, 29, 24, 30, 49 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 14)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 0.5;

                    int[] idx = new int[] { 10, 18, 26, 33, 11, 19, 27, 34, 12, 20, 28, 13, 21, 2, 3, 4, 5 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 15)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 0.5;

                    int[] idx = new int[] { 9, 17, 25, 32, 38, 10, 18, 26, 33, 39, 34, 40, 28, 35, 13, 21, 29, 14, 22, 1, 2, 5, 6 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 16)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 1;

                    int[] idx = new int[] { 10, 18, 26, 27, 40, 45, 28, 35, 29, 2, 50 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 17)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 1;

                    int[] idx = new int[] { 10, 18, 26, 39, 44, 27, 34, 2, 49 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 18)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 0.5;

                    int[] idx = new int[] { 26, 33, 19, 27, 34, 20, 28 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 19)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = -1;

                    int[] idx = new int[] { 10, 18, 26, 33, 39, 44, 2, 49 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 20)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 1;

                    int[] idx = new int[] { 10, 18, 26, 27, 28, 13, 21, 2, 5 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 21)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 0.5;

                    int[] idx = new int[] { 11, 19, 27, 28, 29, 3, 30 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 22)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 1;

                    int[] idx = new int[] { 25, 26, 11, 19, 27, 12, 20, 28, 29, 3, 4, 30, 24 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 23)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 1;

                    int[] idx = new int[] { 25, 43, 26, 44, 11, 19, 40, 45, 35, 29, 14, 22, 3, 6, 7, 15, 23, 30, 36, 41, 46, 47, 48, 49, 50, 51, 24, 42 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 24)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = -1;

                    int[] idx = new int[] { 9, 17, 25, 32, 38, 43, 10, 18, 26, 33, 1, 2, 48, 0, 8, 16, 24, 31, 37, 42, 47 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 25)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = -1;

                    int[] idx = new int[] { 9, 17, 25, 32, 38, 43, 10, 18, 26, 33, 11, 19, 27, 12, 20, 1, 2, 3, 4, 48, 0, 8, 16, 24, 31, 37, 42, 47 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 26)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 1;

                    int[] idx = new int[] { 39, 44, 49 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Selected Cascade size (n_size: " +
                    cas_size.ToString() + ") have not contain case: " + caseNo.ToString() + ". Select between 0 and 26.");
            }
            else if (cas_size == 4)
            {
                if (caseNo == 1)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 4;

                    int[] idx = new int[] { 30, 31 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 2)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 4;

                    int[] idx = new int[] { 38, 31 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 3)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 4;

                    int[] idx = new int[] { 38 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);

                }
                else AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Selected Cascade size (n_size: " +
                    cas_size.ToString() + ") have not contain case: " + caseNo.ToString() + ". Select between 0 and 3.");
            }
            else if (cas_size == 5)
            {
                if (caseNo == 1)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 4;

                    int[] idx = new int[] { 42, 43 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);

                }
                else if (caseNo == 2)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 4;

                    int[] idx = new int[] { 50, 43 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 3)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 4;

                    int[] idx = new int[] { 50, 43, 35 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 4)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 4;

                    int[] idx = new int[] { 33, 50, 43, 35 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 5)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 4;

                    int[] idx = new int[] { 33, 34, 35 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else if (caseNo == 6)
                {
                    Brep sphbrep = sph.ToBrep();
                    double translatedist = 4;

                    int[] idx = new int[] { 50 };
                    for (int i = 0; i < idx.Length; i++)
                    {
                        net_extent[idx[i]] = TranslatePt(sphbrep, net_extent[idx[i]], translatedist);
                    }

                    netpts = NetPtsOut(net_extent, n);
                }
                else AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Selected Cascade size (n_size: " +
                    cas_size.ToString() + ") have not contain case: " + caseNo.ToString() + ". Select between 0 and 6.");
            }
            else AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Selected Cascade size (n_size: " +
                cas_size.ToString() + ") have not contain case: " + caseNo.ToString() + ". Select default case: 0.");
            ///
        }

        private void SetExtNet(int cas_n, List<Point3d> net_extent,out List<Line> extNetLines)
        {
            // Assigning the input pts with respect to their cascade size
            int nn = cas_n + 5;
            Point3d[,] d_ij = new Point3d[nn, nn];
            for (int i = 0, ii = 0, mm = 0; i < nn; i++)
            {
                for (int j = 0; j < nn - mm; j++)
                {
                    d_ij[j, i] = net_extent[ii];
                    ii++;
                }
                for (int j = nn - mm; j < nn; j++)
                {
                    d_ij[j, i] = Point3d.Unset;
                }
                if (i > 1 && i < nn - 3) mm++;
            }

            // Constructing the net edges
            extNetLines = new List<Line>();
            int m = d_ij.GetLength(0), mm0 = 1, mm1 = 0;
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < m - mm1; j++)// Vertical lines
                {
                    if (i == m - 1) break;
                    else if (i == 0 || i == 1) extNetLines.Add(new Line(d_ij[j, i], d_ij[j, i + 1]));
                    else if (i == m - 3 || i == m - 2) extNetLines.Add(new Line(d_ij[j, i], d_ij[j, i + 1]));
                    else
                    {
                        if (j < 3) extNetLines.Add(new Line(d_ij[j, i], d_ij[j, i + 1]));
                        else if (j > (m - mm1) - 4) extNetLines.Add(new Line(d_ij[j, i], d_ij[j - 1, i + 1]));
                        else
                        {
                            extNetLines.Add(new Line(d_ij[j, i], d_ij[j - 1, i + 1]));
                            extNetLines.Add(new Line(d_ij[j, i], d_ij[j, i + 1]));
                        }
                    }
                }
                for (int j = 0; j < m - mm0; j++)// Horizontal lines
                {
                    extNetLines.Add(new Line(d_ij[j, i], d_ij[j + 1, i]));
                }
                if (i > 2 && i < m - 3) mm0++;
                if (i > 1 && i < m - 3) mm1++;
            }
        }

        private List<Point3d> NetPtsOut(List<Point3d> net_extent,int n)
        {
            List<Point3d> netpts = new List<Point3d>();
            for (int i = (n + 2) + 1, ii = 0, iii = 0; i < net_extent.Count - 6; i = i + 2)
            {
                for (int j = 0; j < (n) - iii; j++)
                {
                    netpts.Add(net_extent[i]);
                    i++;
                }
                if (ii > 0 && ii < (n - 2)) iii++;
                ii++;
            }
            return netpts;
        }

        private Point3d SetQMiddle(Point3d p0, Point3d p1, Point3d p2, Point3d p3)
        {
            return (p0 + p1 + p2 + p3) / 4;
        }
        private Point3d SetEMiddle(Point3d p0, Point3d p1)
        {
            return (p0 + p1) / 2;
        }

        private static Point3d TranslatePt(Brep brep, Point3d pt, double translatedist)
        {
            brep.ClosestPoint(pt, out Point3d closestpoint, out ComponentIndex ci, out double s, out double t, 0, out Vector3d normal);
            normal.Unitize();
            return new Point3d(pt + (normal * translatedist));
        }

        private void SetExtendedNCasSurf(int cas_n, List<Point3d> net_extent, out List<Point3d> tensorpoints, out List<NurbsSurface> surfs, out StringBuilder bvfile)
        {
            int n = (cas_n + 3) + 2;

            #region extended surface tensor border construction
            tensorpoints = new List<Point3d>();
            for (int ii = 0; ii < 2; ii++)
            {
                if (ii != 0)
                {
                    for (int i = (ii * n); i < (ii * n) + (n-1); i++)
                    {
                        tensorpoints.Add(SetEMiddle(net_extent[i], net_extent[i + 1]));
                        if (i != (ii * n) + (n-2)) tensorpoints.Add(net_extent[i + 1]);
                    }
                }
                for (int i = (ii * n); i < (ii * n) + (n-1); i++)
                {
                    if (i != (ii * n)) tensorpoints.Add(SetEMiddle(net_extent[i], net_extent[i + n]));
                    tensorpoints.Add(SetQMiddle(net_extent[i], net_extent[i + 1], net_extent[i + n], net_extent[i + (n+1)]));
                }
            }


            List<List<int>> idx = new List<List<int>>();
            for (int i = 0; i < cas_n; i++)
            {
                idx.Add(new List<int>());

                if (i == 0) idx[idx.Count - 1].Add((2 * n));
                else
                {
                    idx[idx.Count - 1].Add(idx[idx.Count - 2][idx[idx.Count - 2].Count - 1] + (cas_n + 3 - i));
                    idx[idx.Count - 1].Add(idx[idx.Count - 1][idx[idx.Count - 1].Count - 1] + 3);
                }
            }

            for (int i = 0; i < cas_n; i++)
            {
                for (int ii = 0; ii < idx[i].Count; ii++)
                {
                    tensorpoints.Add(SetEMiddle(net_extent[idx[i][ii]], net_extent[idx[i][ii] + 1]));
                    tensorpoints.Add(net_extent[idx[i][ii] + 1]);
                    tensorpoints.Add(SetEMiddle(net_extent[idx[i][ii] + 1], net_extent[idx[i][ii] + 2]));

                    tensorpoints.Add(SetQMiddle(net_extent[idx[i][ii]], net_extent[idx[i][ii] + 1], net_extent[idx[i][ii] + ((n) - i)], net_extent[idx[i][ii] + ((n + 1) - i)]));
                    tensorpoints.Add(SetEMiddle(net_extent[idx[i][ii] + 1], net_extent[idx[i][ii] + ((n + 1) - i)]));
                    tensorpoints.Add(SetQMiddle(net_extent[idx[i][ii] + 1], net_extent[idx[i][ii] + 2], net_extent[idx[i][ii] + ((n + 1) - i)], net_extent[idx[i][ii] + ((n + 2) - i)]));
                }
            }

            idx.Add(new List<int>());
            idx[idx.Count - 1].Add(idx[idx.Count - 2][idx[idx.Count - 2].Count - 1] + 3);
            idx[idx.Count - 1].Add(idx[idx.Count - 1][idx[idx.Count - 1].Count - 1] + 3);
            idx[idx.Count - 1].Add(idx[idx.Count - 1][idx[idx.Count - 1].Count - 1] + 2);
            for (int i = idx.Count - 1; i < (idx.Count - 1) + 1; i++)
            {
                for (int ii = 0; ii < idx[i].Count; ii++)
                {
                    tensorpoints.Add(SetEMiddle(net_extent[idx[i][ii]], net_extent[idx[i][ii] + 1]));
                    tensorpoints.Add(net_extent[idx[i][ii] + 1]);
                    tensorpoints.Add(SetEMiddle(net_extent[idx[i][ii] + 1], net_extent[idx[i][ii] + 2]));

                    tensorpoints.Add(SetQMiddle(net_extent[idx[i][ii]], net_extent[idx[i][ii] + 1], net_extent[idx[i][ii] + 5], net_extent[idx[i][ii] + 6]));
                    tensorpoints.Add(SetEMiddle(net_extent[idx[i][ii] + 1], net_extent[idx[i][ii] + 6]));
                    tensorpoints.Add(SetQMiddle(net_extent[idx[i][ii] + 1], net_extent[idx[i][ii] + 2], net_extent[idx[i][ii] + 6], net_extent[idx[i][ii] + 7]));
                    if (ii == 1) tensorpoints.Add(SetEMiddle(net_extent[idx[i][ii] + 2], net_extent[idx[i][ii] + 7]));
                }
            }

            idx.Add(new List<int>());
            idx[idx.Count - 1].Add(idx[idx.Count - 2][idx[idx.Count - 2].Count - 1] + 3);
            idx[idx.Count - 1].Add(idx[idx.Count - 1][idx[idx.Count - 1].Count - 1] + 2);
            for (int i = idx.Count - 1; i < (idx.Count - 1) + 1; i++)
            {
                for (int ii = 0; ii < idx[i].Count; ii++)
                {
                    tensorpoints.Add(SetEMiddle(net_extent[idx[i][ii]], net_extent[idx[i][ii] + 1]));
                    tensorpoints.Add(net_extent[idx[i][ii] + 1]);
                    tensorpoints.Add(SetEMiddle(net_extent[idx[i][ii] + 1], net_extent[idx[i][ii] + 2]));
                    if (ii == 0) tensorpoints.Add(net_extent[idx[i][ii] + 2]);

                    tensorpoints.Add(SetQMiddle(net_extent[idx[i][ii]], net_extent[idx[i][ii] + 1], net_extent[idx[i][ii] + 5], net_extent[idx[i][ii] + 6]));
                    tensorpoints.Add(SetEMiddle(net_extent[idx[i][ii] + 1], net_extent[idx[i][ii] + 6]));
                    tensorpoints.Add(SetQMiddle(net_extent[idx[i][ii] + 1], net_extent[idx[i][ii] + 2], net_extent[idx[i][ii] + 6], net_extent[idx[i][ii] + 7]));
                    if (ii == 0) tensorpoints.Add(SetEMiddle(net_extent[idx[i][ii] + 2], net_extent[idx[i][ii] + 7]));
                }
            }

            #endregion


            #region writing bv file and surface construction

            bvfile = new StringBuilder();
            bvfile.AppendLine("Group 1 extend");
            surfs = new List<NurbsSurface>();
            int nn = (cas_n + 3);
            int nn2 = 2 * nn + 1;
            for (int ii = 0; ii < nn; ii++)
            {
                int i = ii * 2;
                surfs.Add(NurbsSurface.CreateFromPoints(new Point3d[] {
                    tensorpoints[i], tensorpoints[i + 1],  tensorpoints[i + 2],
                    tensorpoints[i + nn2], tensorpoints[i + (nn2+1)], tensorpoints[i + (nn2+2)],
                    tensorpoints[i + (2*nn2)], tensorpoints[i + (2*nn2+1)], tensorpoints[i + (2*nn2+2)] }, 3, 3, 2, 2));

                bvfile.AppendLine("5");
                bvfile.AppendLine("2 2");
                bvfile.AppendLine($"{tensorpoints[i].X.ToString()} {tensorpoints[i].Y} {tensorpoints[i].Z}");
                bvfile.AppendLine($"{tensorpoints[i + 1].X} {tensorpoints[i + 1].Y} {tensorpoints[i + 1].Z}");
                bvfile.AppendLine($"{tensorpoints[i + 2].X} {tensorpoints[i + 2].Y} {tensorpoints[i + 2].Z}");
                bvfile.AppendLine($"{tensorpoints[i + nn2].X} {tensorpoints[i + nn2].Y} {tensorpoints[i + nn2].Z}");
                bvfile.AppendLine($"{tensorpoints[i + (nn2 + 1)].X} {tensorpoints[i + (nn2 + 1)].Y} {tensorpoints[i + (nn2 + 1)].Z}");
                bvfile.AppendLine($"{tensorpoints[i + (nn2 + 2)].X} {tensorpoints[i + (nn2 + 2)].Y} {tensorpoints[i + (nn2 + 2)].Z}");
                bvfile.AppendLine($"{tensorpoints[i + (2 * nn2)].X} {tensorpoints[i + (2 * nn2)].Y} {tensorpoints[i + (2 * nn2)].Z}");
                bvfile.AppendLine($"{tensorpoints[i + (2 * nn2 + 1)].X} {tensorpoints[i + (2 * nn2 + 1)].Y} {tensorpoints[i + (2 * nn2 + 1)].Z}");
                bvfile.AppendLine($"{tensorpoints[i + (2 * nn2 + 2)].X} {tensorpoints[i + (2 * nn2 + 2)].Y} {tensorpoints[i + (2 * nn2 + 2)].Z}");
            }

            int nn3 = (5 + (cas_n * 2));
            int[][] index = new int[nn3][];
            for (int i = 0; i < nn3; i++)
            {
                if (i == 0)
                {
                    int mm = (2 * nn2);
                    index[i] = new int[] { mm, (mm + 1), (mm + 2), (mm + 2) + nn3, (mm + 2) + nn3 + 1, (mm + 2) + nn3 + 2, (mm + 2) + nn3 + 3, (mm + 2) + nn3 + 4, (mm + 2) + nn3 + 5 }; i++;
                    int m = i - 1;
                    int mmm = nn3 - 1;
                    index[i] = new int[] { index[m][0] + mmm, index[m][1] + mmm, index[m][2] + mmm, index[m][3] + 6, index[m][4] + 6, index[m][5] + 6, index[m][6] + 6, index[m][7] + 6, index[m][8] + 6 };
                }
                else if (i >= nn3 - 4)
                {
                    int m = i - 1;
                    index[i] = new int[] { index[m][0] + 6, index[m][1] + 6, index[m][2] + 6, index[m][3] + 7, index[m][4] + 7, index[m][5] + 7, index[m][6] + 7, index[m][7] + 7, index[m][8] + 7 }; i++;
                    m = i - 1;
                    index[i] = new int[] { index[m][0] + 6, index[m][1] + 6, index[m][2] + 6, index[m][3] + 6, index[m][4] + 6, index[m][5] + 6, index[m][6] + 7, index[m][7] + 7, index[m][8] + 7 }; i++;
                    m = i - 1;
                    index[i] = new int[] { index[m][0] + 2, index[m][1] + 2, index[m][2] + 5, index[m][3] + 2, index[m][4] + 2, index[m][5] + 6, index[m][6] + 2, index[m][7] + 2, index[m][8] + 5 }; i++;
                    m = i - 1;
                    index[i] = new int[] { index[m][0] + 5, index[m][1] + 5, index[m][2] + 2, index[m][3] + 6, index[m][4] + 6, index[m][5] + 2, index[m][6] + 5, index[m][7] + 5, index[m][8] + 2 };
                }
                else
                {
                    int m = i - 1;
                    index[i] = new int[] { index[m][0] + 6, index[m][1] + 6, index[m][2] + 6, index[m][3] + 6, index[m][4] + 6, index[m][5] + 6, index[m][6] + 6, index[m][7] + 6, index[m][8] + 6 };
                }
            }

            for (int i = 0; i < index.Length; i++)
            {
                surfs.Add(NurbsSurface.CreateFromPoints(new Point3d[] {
                    tensorpoints[index[i][0]], tensorpoints[index[i][1]], tensorpoints[index[i][2]],
                    tensorpoints[index[i][3]], tensorpoints[index[i][4]], tensorpoints[index[i][5]],
                    tensorpoints[index[i][6]], tensorpoints[index[i][7]], tensorpoints[index[i][8]] }, 3, 3, 2, 2));

                bvfile.AppendLine("5");
                bvfile.AppendLine("2 2");
                bvfile.AppendLine($"{tensorpoints[index[i][0]].X} {tensorpoints[index[i][0]].Y} {tensorpoints[index[i][0]].Z}");
                bvfile.AppendLine($"{tensorpoints[index[i][1]].X} {tensorpoints[index[i][1]].Y} {tensorpoints[index[i][1]].Z}");
                bvfile.AppendLine($"{tensorpoints[index[i][2]].X} {tensorpoints[index[i][2]].Y} {tensorpoints[index[i][2]].Z}");
                bvfile.AppendLine($"{tensorpoints[index[i][3]].X} {tensorpoints[index[i][3]].Y} {tensorpoints[index[i][3]].Z}");
                bvfile.AppendLine($"{tensorpoints[index[i][4]].X} {tensorpoints[index[i][4]].Y} {tensorpoints[index[i][4]].Z}");
                bvfile.AppendLine($"{tensorpoints[index[i][5]].X} {tensorpoints[index[i][5]].Y} {tensorpoints[index[i][5]].Z}");
                bvfile.AppendLine($"{tensorpoints[index[i][6]].X} {tensorpoints[index[i][6]].Y} {tensorpoints[index[i][6]].Z}");
                bvfile.AppendLine($"{tensorpoints[index[i][7]].X} {tensorpoints[index[i][7]].Y} {tensorpoints[index[i][7]].Z}");
                bvfile.AppendLine($"{tensorpoints[index[i][8]].X} {tensorpoints[index[i][8]].Y} {tensorpoints[index[i][8]].Z}");
            }
            #endregion
        }

        public override void AddRuntimeMessage(GH_RuntimeMessageLevel level, string text)
        {
            base.AddRuntimeMessage(level, text);
        }
        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("1c4e5609-e46b-4c98-8d32-60ce77b1f10b"); }
        }
    }
}