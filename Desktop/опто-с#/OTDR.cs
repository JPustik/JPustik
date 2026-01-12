
using System;
using System.Linq;

public static class OTDR
{
    public record Params(
        double P0 = 0.005,
        double tu = 100e-9,
        double Vd = 1.5e5,
        double SMR = 10,
        double K = 1e6,
        double R1 = 0.04,
        double R2 = 0.00003,
        double R3 = 0.04,
        double kps = 0.6,
        bool welded = true,
        double L1 = 1.5,
        double L2 = 2.0,
        double AL1 = 0.46,
        double AL2 = 0.9,
        double sko = 0.008,
        int N = 20000
    );

    static Random rnd = new();

    static double Gaussian(double mean, double std)
    {
        var u1 = 1.0 - rnd.NextDouble();
        var u2 = 1.0 - rnd.NextDouble();
        var randStd = Math.Sqrt(-2.0 * Math.Log(u1)) *
                      Math.Sin(2.0 * Math.PI * u2);
        return mean + std * randStd;
    }

    public static (double[] distance, double[] logF) Calculate(Params p)
    {
        double AS1 = p.AL1 * 4e-4;
        double AS2 = p.AL2 * 4e-4;

        double tt1 = 2 * p.L1 / p.Vd - p.tu / 2;
        double tt2 = tt1 + 2 * p.L2 / p.Vd;

        double dt = p.tu / 10;
        double tmax = tt2 + 5 * p.tu;
        int n = (int)(tmax / dt);

        double[] t = Enumerable.Range(0, n).Select(i => i * dt).ToArray();
        double[] F = new double[n];

        for (int i = 0; i < n; i++)
        {
            double ti = t[i];
            double Fi;
            if (ti <= p.tu)
                Fi = -p.P0 * (p.R1 + AS1 * p.Vd * p.tu * ti);
            else if (ti <= tt1)
                Fi = -p.P0 * (p.R1 * Math.Exp(-3 * (ti - p.tu) / p.tu));
            else
                Fi = -p.P0 * p.R1 * Math.Exp(-3 * (ti - p.tu) / p.tu);

            F[i] = p.K * p.SMR * Fi;
        }

        double noiseStd = p.sko / Math.Sqrt(p.N);
        double[] logF = new double[n];
        double[] dist = new double[n];

        for (int i = 0; i < n; i++)
        {
            double noisy = F[i] + Gaussian(0, noiseStd);
            logF[i] = noisy < 0 ? 5 * Math.Log10(Math.Max(-noisy, 1e-15) / 5) : -60;
            dist[i] = p.Vd * t[i] / 2;
        }

        return (dist, logF);
    }
}
