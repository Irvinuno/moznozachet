import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import spherical_jn, spherical_yn
c = 3e8
class VariantLoader:
    def __init__(self, filename, variant_number):
        self.filename = filename
        self.variant_number = variant_number
    def load(self):
        with open(self.filename, "r", encoding="utf-8") as f:
            data = json.load(f)["data"]
        for item in data:
            v = item["variant"]
            if v["number"] == self.variant_number:
                return {
                    "D": float(v["D"]),
                    "fmin": float(v["fmin"]),
                    "fmax": float(v["fmax"]),
                }
        raise ValueError("Вариант не найден")
class RCSCalculator:
    def __init__(self, radius):
        self.r = radius
    def hankel(self, n, x):
        return spherical_jn(n, x) + 1j * spherical_yn(n, x)
    def calc_sigma(self, f, n_max=30):
        lam = c / f
        k = 2 * np.pi / lam
        kr = k * self.r
        s = 0
        for n in range(1, n_max + 1):
            jn = spherical_jn(n, kr)
            jn_1 = spherical_jn(n - 1, kr)
            yn = spherical_yn(n, kr)
            yn_1 = spherical_yn(n - 1, kr)
            hn = jn + 1j * yn
            hn_1 = jn_1 + 1j * yn_1
            a_n = jn / hn
            b_n = (kr * jn_1 - n * jn) / (kr * hn_1 - n * hn)
            s += ((-1) ** n) * (n + 0.5) * (b_n - a_n)
        sigma = (lam ** 2 / np.pi) * abs(s) ** 2
        return sigma.real
class ResultWriter:
    def __init__(self, filename):
        self.filename = filename
    def write(self, freq, sigma):
        with open(self.filename, "w", encoding="utf-8") as f:
            for f_i, s_i in zip(freq, sigma):
                f.write(f"{f_i:.6e}    {s_i:.6e}\n")
class Plotter:
    @staticmethod
    def plot(freq, sigma):
        plt.figure(figsize=(8, 5))
        plt.plot(freq / 1e9, sigma)
        plt.xlabel("Частота, ГГц")
        plt.ylabel("ЭПР, м²")
        plt.grid(True)
        plt.tight_layout()
        plt.show()
def main():
    loader = VariantLoader("task_rcs_01.json", 16)
    params = loader.load()
    r = params["D"] / 2
    fmin, fmax = params["fmin"], params["fmax"]
    freq = np.linspace(fmin, fmax, 500)
    calculator = RCSCalculator(r)
    sigma = np.array([calculator.calc_sigma(f) for f in freq])
    writer = ResultWriter("rcs_variant_16.txt")
    writer.write(freq, sigma)
    Plotter.plot(freq, sigma)
if __name__ == "__main__":
    main()
