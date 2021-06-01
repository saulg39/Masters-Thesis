# Masters-Thesis

I_beam.py					defines I-Beam node coordinates and connections(strips)
I_beam_res.py					defines I-beam node coordinates and connections(strips) and reidual stresses
Moment_elastic_crit.py				Calculates elastic buckling moment for working out slenderness(not used in end)
RHS.py						defines RSH node coordinates and connections(strips)			
channel.py					defines channel node coordinates and connections(strips)
finite_strip_Compression_new.py			FSM for compression
finite_strip_compression_curve_new.py		Draws curve from compression FSM
finite_strip_curvature_curve.py			Draw curvature against length curve for moment FSM
finite_strip_elastic_stress.py			Calculates elastic buckling stress for working out slenderness
finite_strip_moment.py				FSM for constent bending moment
finite_strip_moment_curve.py			Draws curve from moment FSM
finite_strip_moment_curve_res.py		Draws curve from moment FSM taking in to account residual stresses
finite_strip_moment_nweb.py			varies number of stripes in web so optimal amount can be found
finite_strip_moment_output.py			writes outputes in to excel sheet
finite_strip_moment_res.py			FSM for moment taking in to account residual stresses
get_data.py					retrieves data from excel
graph_code.py					Code for data retrival (very complex)
graph_code_2.py					Code for data retrival (very complex)
graph_moment.py					Calculate theoretical moment from graph
mid_axes.py					Calculate N.A. of shape (not used due to symmetry of all shapes)
plastic_moment.py				Calculate plastic moment used for the data retreval
stress_from_strain.py				returns stress and Et from strain and material properties

Look in Data file to find beam_column_data.xlsx which stores all the used data and outcomes from CSM and Eurocode (used to help make graphs in Latex so is a bit messy)

