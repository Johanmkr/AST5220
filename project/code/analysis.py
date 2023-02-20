import plotly.io as pio
pio.kaleido.scope.mathjax = None
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pandas as pd
pd.options.plotting.backend = "plotly"
from xhtml2pdf import pisa
from IPython import embed

temp_output_path = "../output/"

def convert_html_to_pdf(source_html, output_filename):
    result_file = open(temp_output_path + output_filename, "w+b")

    pisa_status = pisa.CreatePDF(
        source_html,
        dest=result_file
    )

    result_file.close()

    return pisa_status.err

# Using numpy
supernova_data = np.loadtxt("data/supernovadata.txt")
z = supernova_data[:,0]
x_of_a = -np.log(1+z)
d_L = supernova_data[:,1]
error = supernova_data[:,2]

# Using pandas
supernova_frame =pd.read_fwf("data/supernovadata.txt").drop(columns=["#"])

header_names = ["x", "eta", "Hp", "dHp", "OmegaB", "OmegaCDM", "OmegaLambda", "OmegaR", "OmegaNu", "OmegaK", "d_L", "dummy"]

cosmology = pd.read_csv("data/backgroundcosmology-2_0.txt", delimiter=" ", names=header_names)

# embed()
# Create figure of luminosity distance

fig_lumdist = go.Figure()
fig_lumdist.add_trace(go.Scatter(x=cosmology["x"], y=cosmology["d_L"]/1e9,
                              name = "Luminosity distance"))
#   add error
fig_lumdist.add_trace(go.Scatter(x=x_of_a, y=d_L,
                                 name="Error",
                                 mode="markers",
                                 error_y=dict(
                                    array=error,
                                    color="red"
                                 )))

fig_lumdist.update_layout(
    title="Supernova data",
    showlegend=True,
    xaxis_title=r"$x$",
    # xaxis_range=[-2,0],
    # yaxis_range=[0,np.max(d_L)],
    yaxis_title=r"$d_L$",
    legend_title="Legend",
    font=dict(
        family="Times New Roman",
        size=18,
        color="black"
    ))
fig_lumdist.update_yaxes(type="log")
fig_lumdist.show()

# fig_lumdist.write_image(temp_output_path+"lumdist.pdf")

# fig_lumdist.write_html(temp_output_path+"lumdist.html")

# convert_html_to_pdf(temp_output_path + "lumdist.html", "lumdist.pdf")

# fig_lumdist=make_figure(df, pa)

# pl.io.write_image(fig_lumdist, temp_output_path+"lumdist.pdf", format="pdf")







if __name__=="__MAIN__":
    pass