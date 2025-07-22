from dash import Dash, html, dcc, Input, Output, State, callback, ctx
import plotly.graph_objs as go
from anytree import Node, RenderTree
import base64

# Create the Dash app
app = Dash(__name__)

polyline = [ (0, 0),
(2, -1),
(4, 1),
(3, 2),
(2, 1),
(4, -1),
(6, 0) ]

simplification = [0, 6]
epsilon = 2 
distance_measure = 'E'


data = {
  (0, 0, 0): (0, 0, 0, 0),

  (1, 1, 0): (0, 0, 0, 0.1056),
  (1, 1, 1): (0, 0, 0, 0),
  (1, 2, 1): (0, 0, 0, 0.2928),
  (1, 2, 2): (0, 0, 0, 0),
  (1, 2, 3): (0, 0, 0, 0),
  (1, 2, 4): (0, 0, 0, 0),
  (1, 2, 5): (0, 0, 0, 0),
  (1, 3, 1): (0, 0, 0, 0.5),
  (1, 3, 2): (0, 0, 0, 0),
  (1, 3, 3): (0, 0, 0, 0),
  (1, 3, 4): (0, 0, 0, 0),
  (1, 4, 0): (0, 0, 0, 0.2),
  (1, 4, 1): (0, 0, 0, 0),
  (1, 4, 2): (0, 0, 0, 0),
  (1, 4, 3): (0, 0, 0, 0),
  (1, 4, 4): (0, 0, 0, 0),
  (1, 5, 0): (0, 0, 0, 1),
  (1, 5, 1): (0, 0, 0, 0),
  (1, 5, 2): (0, 0, 0, 0),
  (1, 6, 5): (0, 0, 0, 0.1056)
}

def read_file(file):
    global data, simplification, epsilon, polyline, distance_measure
    lines = file.strip().split('\n')
    first_line = lines[0].strip().split(' ')
    if len(first_line) != 4:
        raise ValueError("First line must contain exactly 4 values: Distance, epsilon, n, d")

    # Parse Distance (either E, M, C, IE, or an integer > 2)
    D = first_line[0]
    if D not in ['E', 'M', 'C', 'IE']:
        try:
            D = int(D)
            if D <= 2:
                raise ValueError("Distance must be either E, M, C, IE, or an integer > 2")
        except ValueError:
            raise ValueError("Distance must be either E, M, C, IE, or an integer > 2")
    distance_measure = D

    # Parse epsilon (float > 0)
    try:
        epsilon = float(first_line[1])
        if epsilon < 0:
            raise ValueError("epsilon must be a float >= 0")
    except ValueError:
        raise ValueError("epsilon must be a float >= 0")

    # Parse n and d
    try:
        n = int(first_line[2])
        d = int(first_line[3])
        if n <= 0 or d <= 0:
            raise ValueError("n and d must be positive integers")
    except ValueError:
        raise ValueError("n and d must be integers")

    if d != 2:
        raise ValueError("Dimensions other than 2 are not supported")

    # Parse the next n lines (each with d floats)
    polyline = []
    for i in range(n):
        line = lines[i + 1].strip().split(' ')
        if len(line) != d:
            raise ValueError(f"Expected {d} floats, but found {len(line)}")
        try:
            point = tuple(map(float, line))
            polyline.append(point)
        except ValueError:
            raise ValueError("All values in the polyine must be floats")


    # Parse the remaining lines simplification data
    k = 0
    data = {}
    for line in lines[1 + n: ]:
        line = line.strip().split(' ')
        if len(line) != 6:
            raise ValueError("Each additional line must contain 5 integers and 1 float or integer in the implicit case")
        try:
            integers = list(map(int, line[:5]))
            if D == 'IE':
                t = int(line[5])
                if t < 0:
                    raise ValueError("the point reference must be non negative")
            else:
                t = float(line[5])
                if t < 0 or t > 1:
                    raise ValueError("t must be between 0 and 1")
            data[integers[0], integers[1], integers[2]] = (max(integers[0] - 1, 0), integers[3], integers[4], t)
            k = max(k, integers[0])
        except ValueError:
            raise ValueError("Invalid format in additional lines")

    i = n - 1 
    j = n - 2 
    simpl_values = [i]
    while k > 0:
        k, i, j, _ = data[k, i, j]
        simpl_values.append(i)
    simplification = list(reversed(simpl_values))



tree_nodes = {}
def make_tree():
    global tree_nodes

    tree_nodes = {}
    for d, v in data.items():
        name = str(d[0]) + "," + str(d[1]) + "," + str(d[2])
        if d[0] == 0:
            if d[2] == 0:
                tree_nodes[d] = Node(name, data=dict(k=d[0], i=d[1], j=d[2]))
            else:
                tree_nodes[d] = Node(name, parent=tree_nodes[0, 0, d[2] - 1], data=dict(k=d[0], i=d[1], j=d[2]))
        else:
            tree_nodes[d] = Node(name, parent=tree_nodes[v[0], v[1], v[2]], data=dict(k=d[0], i=d[1], j=d[2]))

make_tree()

lookup_point_name_dict = {}


distance_shape_index = [-1] * len(polyline)
distance_shown = [False] * len(polyline)
distance_shape_layouts = []

def reset_distances():
    global distance_shape_layouts, distance_shown, distance_shape_index, distance_measure
    distance_shape_index = [-1] * len(polyline)
    distance_shown = [False] * len(polyline)
    distance_shape_layouts = []
    for v in polyline: 
        if distance_measure == 'E' or distance_measure == 'IE':
            distance_shape_layouts.append({
                'type': "circle",
                'xref': "x",
                'yref': "y",
                'x0': v[0] - epsilon,
                'y0': v[1] - epsilon,
                'x1': v[0] + epsilon,
                'y1': v[1] + epsilon,
                'fillcolor': "lightgreen", 
                'opacity':0, 
                'line':dict(color="blue")
            })
        elif distance_measure == 'C':
            distance_shape_layouts.append({
                'type': "path",
                'path': f'M {v[0]-epsilon}, {v[1]-epsilon} L {v[0]+epsilon}, {v[1]-epsilon} L {v[0]+epsilon}, {v[1]+epsilon} L {v[0]-epsilon}, {v[1]+epsilon} Z',
                'xref': "x",
                'yref': "y",
                'fillcolor': "lightgreen", 
                'opacity':0, 
                'line':dict(color="blue")
            })
        elif distance_measure == 'M':
            distance_shape_layouts.append({
                'type': "path",
                'path': f'M {v[0]-epsilon}, {v[1]} L {v[0]}, {v[1]-epsilon} L {v[0]+epsilon}, {v[1]} L {v[0]}, {v[1]+epsilon} Z',
                'xref': "x",
                'yref': "y",
                'fillcolor': "lightgreen", 
                'opacity':0, 
                'line':dict(color="blue")
            })



reset_distances()

def draw_tree(node, layer = 0, x = 0, opacity=1):
    siblings = [c for c in node.children if c.data['k'] == node.data['k']]
    children = [c for c in node.children if c.data['k'] != node.data['k']]

    if children:
        total_width = -1
    else:
        total_width = 0 

    total_data = []
    child_pos = []
    for c in children:
        width, data, c_x, c_y, name = draw_tree(c, layer + 1, x + total_width + 1)
        total_data.extend(data)
        total_width += width + 1
        child_pos.append((c_x, c_y, name))


    pos_x = x + total_width / 2
    pos_y = layer

    for _x, _y, cname in child_pos:
        total_data.append(go.Scatter(x=[pos_x, _x],y=[pos_y, _y], text=["l" + cname], mode='lines', opacity=opacity, line=dict(color="brown", width=2)))

    for s in siblings:
        width, data, c_x, c_y, name = draw_tree(s, layer, x + total_width + 1)
        total_data.extend(data)
        total_width += width + 1

        total_data.append(go.Scatter(x=[pos_x, c_x], y=[pos_y, c_y], text=["l"+name], mode='lines', opacity=opacity, line=dict(color="brown", width=2)))
        child_pos.append((c_x, c_y, name))

    for x, y, child_name in child_pos:
        total_data.append(go.Scatter(x=[x],y=[y], text=[child_name], mode='markers+text', opacity=opacity, marker=dict(color="green", size=10), textposition="top center"))

    if layer == 0:
        total_data.append(go.Scatter(x=[pos_x],y=[pos_y], text=[node.name], opacity=opacity, mode='markers+text', marker=dict(color="green", size=10), textposition="top center"))

    return total_width, total_data, pos_x, pos_y, node.name

polyfig = None
tree_fig = None

def create_tree_fig(name=None):
    global tree_fig
    if name:
        _, tree_traces, _, _, _ = draw_tree(tree_nodes[0,0,0],0,0,0.2)
        [node] = [t for t in tree_traces if t["text"] == name]
    else:
        _, tree_traces, _, _, _ = draw_tree(tree_nodes[0,0,0],0,0,1)
    tree_fig = go.Figure(data=tree_traces, layout=go.Layout(showlegend=False, title="Simplification Tree"))

    tree_fig.update_xaxes(
        showticklabels=False,
        showgrid=False,
        zeroline=False,
        visible=False
    )

    tree_fig.update_yaxes(
        title="k = Amount of line segments",
        title_font=dict(size=18, color='red'),
        tickfont=dict(size=14, color='blue'),
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray',
        zeroline=True,
        zerolinewidth=2,
        zerolinecolor='black',
        dtick=1,
        tickmode="linear"
    )

    return tree_fig


tree_fig = create_tree_fig()

# Layout of the app
app.layout = html.Div([
    dcc.Graph(id='tree-plot', figure=tree_fig, config=dict(doubleClick="reset+autosize")),
    html.Button('Show Simplification', id='back-button', n_clicks=0),
    dcc.Upload(id='upload-file', children=html.Button('Select Visualization File'), multiple=False),
    dcc.Graph(id='polyline-plot', style=dict(width="100vw", height="100vw"))
])

def global_polyline():
    global polyfig
    fig = go.Figure()
    x = [i for i, _ in polyline]
    y = [i for _, i in polyline]
    fig.add_trace(go.Scatter(x=x, y=y, mode='lines', name='Polyline'))
    polyfig = fig
    return fig

def add_simplification(fig):
    simpl = [polyline[i] for i in simplification]
    x_simpl = [i for i, _ in simpl]
    y_simpl = [i for _, i in simpl]
    fig.add_trace(go.Scatter(x=x_simpl, y=y_simpl, mode='lines', name='Simplification', hoverinfo="skip"))
    return fig

def compute_inner_point(j, t):
    global polyline 
    [ux, uy] = polyline[j]
    [vx, vy] = polyline[j+1]
    return [(1 - t) * ux + t * vx, (1 - t) * uy + t * vy]

currfile = None

@app.callback(
    Output('polyline-plot', 'figure'),
    Output('tree-plot', 'clickData'),
    Output('polyline-plot', 'clickData'),
    Output('tree-plot', 'figure'),
    Input('tree-plot', 'clickData'),
    Input('back-button', 'n_clicks'),
    Input('polyline-plot', 'clickData'),
    Input('upload-file', 'filename'),
    Input('upload-file', 'contents'),
)
def update_polyline(click_data, n_clicks, polyline_click, filename, content):
    global epsilon, polyfig, distance_shown, distance_shape_layouts, distance_shape_index, tree_fig, data, polyline, lookup_point_name_dict, distance_measure

    if ctx.triggered_id == 'upload-file':
        _, content_string = content.split(',')
        decoded = base64.b64decode(content_string)
        text_data = decoded.decode("utf-8")

        read_file(text_data)
        make_tree()
        
        fig1 = add_simplification(global_polyline())
        fig2 = create_tree_fig()
        reset_distances()
        return fig1, None, None, fig2

    if ctx.triggered_id == 'back-button':
        fig = global_polyline()
        reset_distances()
        for a in tree_fig.data:
            a.opacity = 1
        
        return add_simplification(fig), None, None, tree_fig

    if ctx.triggered_id == 'polyline-plot' and polyline_click is not None:
        #if "text" in polyline_click['points'][0]:
        #    index = lookup_point_name_dict[polyline_click['points'][0]['text']]
        #else:
        index = polyline_click['points'][0]['pointIndex'] 
        fig = polyfig

        if distance_shape_index[index] == -1:
            distance_shape_index[index] = len(fig.layout.shapes)
            fig.add_shape(distance_shape_layouts[index])

        distance_shown[index] = not distance_shown[index]
        if distance_shown[index]:
            fig.layout.shapes[distance_shape_index[index]]['opacity'] = 0.25
        else:
            fig.layout.shapes[distance_shape_index[index]]['opacity'] = 0


        return fig, None, None, tree_fig

    if click_data is None or 'points' not in click_data or not click_data['points']:
        fig = global_polyline()
        return add_simplification(fig), None, None, tree_fig
    else:
        reset_distances()

    # Get the name of the clicked node
    node_name = click_data['points'][0]['text']
    for a in tree_fig.data:
        a.opacity = 0.2
    tree_fig.data[click_data['points'][0]["curveNumber"]].opacity = 1
    fig = global_polyline()

    k,i,j = map(int, node_name.split(","))
    i_ = 0 
    j_ = 0 
    t_ = 0 
    t  = 0
    if k:
        k_, i_, j_, t = data[k,i,j]
        _, _, _, t_ = data[k_, i_, j_]

    # subpolyline form t' to t 
    sub_poly_points = []
    if distance_measure != 'IE':
        if t_ == 1:
            sub_poly_points.append(polyline[j_ + 1])
        else:
            sub_poly_points.append(compute_inner_point(j_, t_))
            if j_ + 1 < j:
                sub_poly_points.append(polyline[j_ + 1])

    for x in range(j_ + 2, j):
        sub_poly_points.append(polyline[x])

    if distance_measure != 'IE':
        if t == 0:
            sub_poly_points.append(polyline[j])
        else:
            sub_poly_points.append(polyline[j])
            sub_poly_points.append(compute_inner_point(j, t))

    x = [a[0] for a in sub_poly_points]
    y = [a[1] for a in sub_poly_points]
    if distance_measure != 'IE':
        fig.add_trace(go.Scatter(x=x, y=y, mode='lines+markers', name="Subpolyline P[t'...t]", hoverinfo="skip"))
    else:
        fig.add_trace(go.Scatter(x=x, y=y, mode='lines+markers', name="Subpolyline P[j'+1...j]", hoverinfo="skip"))

    simplification_points = [i]
    tk = k
    ti = i 
    tj = j
    while tk > 0:
        [temp] = [a for a in tree_fig.data if a.text is not None and a.text[0][0] == 'l' and list(map(int, a.text[0][1:].split(","))) == [tk, ti, tj]]
        temp.opacity = 1
        tk, ti, tj, _ = data[tk, ti, tj]
        simplification_points.append(ti)

        [temp] = [a for a in tree_fig.data if a.text is not None and a.text[0][0] != 'l' and list(map(int, a.text[0].split(","))) == [tk, ti, tj]]
        temp.opacity = 1

    x = [polyline[a][0] for a in simplification_points]
    y = [polyline[a][1] for a in simplification_points]
    # tree_fig.data[].opacity = 1

    fig.add_trace(go.Scatter(x=x, y=y, mode='lines+markers', name="Simplification of P[0...i]", hoverinfo="skip"))

        

    point_dict = {}
    new_points = [("i", i), ("i'", i_), ("j", j), ("j'", j_), ("j+1", j+1), ("j'+1", j_+1)]
    if t == 0 or t == 1 and distance_measure != 'IE':
        new_points.append(("t", j + t))
    if t_ == 0 or t_ == 1 and distance_measure != 'IE':
        new_points.append(("t'", j_ + t_))
    if distance_measure == 'IE':
        new_points.append(("r'", t_))
        new_points.append(("r", t))
    for x, y in new_points:
        if y in point_dict:
            point_dict[y] = point_dict[y] + ", " + x
        else: 
            point_dict[y] = x


    lookup_point_name_dict = {}
    for x, y in point_dict.items():
        fig.add_trace(go.Scatter(x=[polyline[x][0]], y=[polyline[x][1]], mode="markers+text", text=[y], opacity=1, marker=dict(color="red", size=10), textposition='top center', name=y, hoverinfo="skip"))
        lookup_point_name_dict[y] = x

    fig.update_layout(title=f"Polyline for {node_name}", xaxis_title="X", yaxis_title="Y")
    return fig, None, None, tree_fig

# Run the app
if __name__ == '__main__':
    app.run(debug=True)
