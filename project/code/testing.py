import plotly.io as pio
pio.kaleido.scope.mathjax = None
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pandas as pd
pd.options.plotting.backend = "plotly"
from xhtml2pdf import pisa
from IPython import embed
import astropy.constants as const
from astropy import units

temp_output_path = "../output/"

header_names = ["x", "eta", "Hp", "dHp", "OmegaB", "OmegaCDM", "OmegaLambda", "OmegaR", "OmegaNu", "OmegaK", "dummy"]

cosmology = pd.read_csv("data/backgroundcosmology.txt", delimiter=" ", names=header_names)

# print(cosmology)


def testing_Omegas():
    xvals = cosmology["x"]
    OmegaRel = cosmology["OmegaR"] + cosmology["OmegaNu"]
    OmegaM = cosmology["OmegaB"] + cosmology["OmegaCDM"]
    omegaplot = go.Figure()
    #   Add OmegaRel
    omegaplot.add_trace(
        go.Scatter(
            x=xvals,
            y=OmegaRel,
            mode="lines",
            name=r"$\mathrm{r}$",
            line=dict(
                color="red"
            )
        )
    )
    #   Add OmegaM
    omegaplot.add_trace(
        go.Scatter(
            x=xvals,
            y=OmegaM,
            mode="lines",
            name=r"$\mathrm{M}$ ",
            line=dict(
                color="orange"
            )
        )
    )
    #   Add OmegaLambda
    omegaplot.add_trace(
        go.Scatter(
            x=xvals,
            y=cosmology["OmegaLambda"],
            mode="lines",
            name=r"$\Lambda$ ",
            line=dict(
                color="green"
            )
        )
    )

    #   Add sum
    omegaplot.add_trace(
        go.Scatter(
            x=xvals,
            y=OmegaM+OmegaRel+cosmology["OmegaLambda"],
            mode="lines",
            name="Sum",
            line=dict(
                color="blue"
            )
        )
    )

    #   Update layout
    omegaplot.update_layout(
        title=r"$\Omega_i(x)$ ",
        # showlegend=True,
        xaxis_title=r"$x$",
        legend_title=r"$\Omega_i$",
        font=dict(
            family="Times New Roman",
            size=18,
            color="black"
        ),
        legend=dict(
            # title=r"$\Omega_i$",
            # orientation="h",
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01
        )
    )

    #   SHOW
    omegaplot.show()

#   Plot of eta(x)Hp(x)/c
#   Needs more work
def eta_of_x_Hp():
    etaHpplot = go.Figure()
    etaHpplot.add_trace(
        go.Scatter(
            x=cosmology["x"],
            y=cosmology["eta"]*cosmology["Hp"] / const.c,
            mode="lines",
            name="some",
            line=dict(
                color="blue"
            )
        )
    )
    etaHpplot.update_layout(
        title="sometitle",
        xaxis_title=r"$x$",
        xaxis_range=[-14,0],
        yaxis_range=[.5,4],
        showlegend=True,
        font=dict(
            family="Times New Roman",
            size=18,
            color="black"
        )
    )
    etaHpplot.show()

def eta_of_x():
    etaplot = go.Figure()
    etaplot.add_trace(
        go.Scatter(
            x=cosmology["x"],
            y=cosmology["eta"]*units.m.to("Mpc"),
            mode="lines",
            name="some",
            line=dict(
                color="blue",
            )
        )
    )
    etaplot.update_yaxes(
        type="log"
    )
    etaplot.update_layout(
        title="sometitle",
        xaxis_title=r"$x$",
        xaxis_range=[-12,0],
        showlegend=True,
        font=dict(
            family="Times New Roman",
            size=18,
            color="black"
        )
    )
    etaplot.show()

testing_Omegas()
eta_of_x_Hp()
eta_of_x() 
# if __name__=="__MAIN__":
#     testing_Omegas()