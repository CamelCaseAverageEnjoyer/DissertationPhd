import plotly as plt
import plotly.express as px
import plotly.graph_objs as go

def plot_all(o):
    fig = go.Figure()
    fig.add_scatter(x=o.x, y=o.y)
    fig.show()
