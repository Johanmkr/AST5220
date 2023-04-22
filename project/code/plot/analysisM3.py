from plot_utils import *

# class Perturbation:
#     def __init__(self, list:DataList)->None:
#         self.DataList = DataList

class Perturbation(Milestone):
    def __init__(self, data:Data)->None:
        super().__init__(data)


    def monopole(self)->None:
        fig, ax = plt.subplots()

        ax.plot(self.x, self.T0)
        
        fig.show()


if __name__=="__main__":
    Pert = Perturbation(Data("perturbations_k0.01.csv"))
    Pert.monopole()