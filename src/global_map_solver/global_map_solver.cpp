// Copyright (C) 2018 by Pierre-Yves Lajoie <lajoie.py@gmail.com>

#include "global_map_solver/global_map_solver.h"
#include <math.h>
#include <sstream>      // std::istringstream

namespace global_map_solver {

const std::string GlobalMapSolver::CONSISTENCY_MATRIX_FILE_NAME = std::string("results/consistency_matrix.clq.mtx");
const std::string GlobalMapSolver::CONSISTENCY_LOOP_CLOSURES_FILE_NAME = std::string("results/consistent_loop_closures.txt");

GlobalMapSolver::GlobalMapSolver(const robot_local_map::RobotLocalMap& robot1_local_map,
                const robot_local_map::RobotLocalMap& robot2_local_map,
                const robot_local_map::RobotMeasurements& interrobot_measurements): 
                pairwise_consistency_(robot1_local_map.getTransforms(), robot2_local_map.getTransforms(), 
                            interrobot_measurements.getTransforms(), interrobot_measurements.getLoopClosures(),
                            robot1_local_map.getTrajectory(), robot2_local_map.getTrajectory(),
                            robot1_local_map.getNbDegreeFreedom()){}

SESync::measurements_t GlobalMapSolver::fillMeasurements(const std::vector<int>& max_clique_data){

    // Preallocate output vector
    SESync::measurements_t measurements;

    for (auto const& t : pairwise_consistency_.getTransformsRobot1().transforms)
    {
        // std::cout << t.second.pose.pose.position.x;
        measurements.push_back(graph_utils::convertTransformToRelativePoseMeasurement(t.second));
    }

    for (auto const& t : pairwise_consistency_.getTransformsRobot2().transforms)
    {
        measurements.push_back(graph_utils::convertTransformToRelativePoseMeasurement(t.second));
    }

    for (auto const& t : pairwise_consistency_.getTransformsInterRobot().transforms)
    {
        measurements.push_back(graph_utils::convertTransformToRelativePoseMeasurement(t.second));
    }

    return measurements;
}

int GlobalMapSolver::solveGlobalMap() {
    // Compute consistency matrix
    Eigen::MatrixXi consistency_matrix = pairwise_consistency_.computeConsistentMeasurementsMatrix();
    // graph_utils::printConsistencyGraph(consistency_matrix, CONSISTENCY_MATRIX_FILE_NAME);

    grpahIOExt gioExt; 
    FMC::CGraphIO gio;

    if(gioExt.readGraphMtx(consistency_matrix)){
        gio.m_vi_Vertices.assign(gioExt.m_vi_Vertices.begin(), gioExt.m_vi_Vertices.end());
        gio.m_vi_Edges.assign(gioExt.m_vi_Edges.begin(), gioExt.m_vi_Edges.end());
        // gio.m_vd_Values.assign(gioExt.m_vd_Values.begin(), gioExt.m_vd_Values.end());
        gio.CalculateVertexDegrees();
    }
    // Compute maximum clique
    int max_clique_size = 0;
    std::vector<int> max_clique_data;
    max_clique_size = FMC::maxClique(gio, max_clique_size, max_clique_data);
    // std::cout << "dwq test: max_clique_data:" << max_clique_data[1] << std::endl; 

    // Print results
    graph_utils::printConsistentLoopClosures(pairwise_consistency_.getLoopClosures(), max_clique_data, CONSISTENCY_LOOP_CLOSURES_FILE_NAME);
    // Clean up
    max_clique_data.clear();

    // Fill measurements
    SESync::measurements_t measurements = fillMeasurements(max_clique_data);
    
    // SE-Sync options
    SESync::SESyncOpts opts;
    opts.verbose = true;
    opts.num_threads = 4;

    /// RUN SE-SYNC! (optimization)
    SESync::SESyncResult results = SESync::SESync(measurements, opts);    

    return max_clique_size;
}

bool grpahIOExt::readGraphMtx( const Eigen::MatrixXi& consistency_matrix, float connStrength){
    // char data_type[LINE_LENGTH] = "pattern";
	// char storage_scheme[LINE_LENGTH] = "symmetric";
    int col=0, row=0, rowIndex=0, colIndex=0;
    int entry_counter = 0, num_of_entries = 0;
    
    // Format edges.
    std::stringstream ss;
    for (int i = 0; i < consistency_matrix.rows(); i++) {
      for (int j = i; j < consistency_matrix.cols(); j++) {
		// std::cout << consistency_matrix(i,j) << ","; 
        if (consistency_matrix(i,j) == 1) {
          ss << i+1 << " " << j+1 << std::endl;
          num_of_entries++;
        }
      }
    }

    std::cout << "num_of_entries: " << num_of_entries << std::endl;
    row = consistency_matrix.rows();
    col = consistency_matrix.cols();
    if(row!=col) 
	{
		cout<<"* WARNING: GraphInputOutput::ReadMatrixMarketAdjacencyGraph()"<<endl;
		cout<<"*\t row!=col. This is not a square matrix. Can't process."<<endl;
		return false;
	}
    
    istringstream in(ss.str());
    istringstream in2; 
    string line = ""; 
    double value; 
    map<int,vector<int> > nodeList;
	map<int,vector<double> > valueList;
	bool b_getValue = true;
	int num_upper_triangular = 0;
    while(!in.eof() && entry_counter<num_of_entries) //there should be (num_of_entries+1) lines in the input file (excluding the comments)
	{
		std::getline(in,line);
		entry_counter++;


		if(line!="")
		{
            
			in2.clear();
			in2.str(line);

			in2 >> rowIndex >> colIndex >> value;
            // std::cout << "row, col, value: "<< rowIndex << " " << colIndex << " " << value << std::endl;
			rowIndex--;
			colIndex--;

			if(rowIndex < 0 || rowIndex >= row)
				cout << "Something wrong rowIndex " << rowIndex << " row " << row << endl;

			if(colIndex < 0 || colIndex >= col)
				cout << "Something wrong colIndex " << colIndex << " col " << col << endl;

			if(rowIndex == colIndex)
			{
				continue;
			}

			// This is to handle directed graphs. If the edge is already present, skip. If not add.
			int exists=0;
			for(int k=0; k<nodeList[rowIndex].size(); k++) {
				if(colIndex == nodeList[rowIndex][k]) {
					exists = 1;
					break;
				}
			}

			if(exists==1) {
				num_upper_triangular++;
			} else {
				if(b_getValue)
				{
					if(value > connStrength)
					{
						nodeList[rowIndex].push_back(colIndex);
						nodeList[colIndex].push_back(rowIndex);
					}
				} 
				else 
				{
					nodeList[rowIndex].push_back(colIndex);
					nodeList[colIndex].push_back(rowIndex);
				}

				if(b_getValue && value > connStrength) 
				{
					valueList[rowIndex].push_back(value);
					valueList[colIndex].push_back(value);
				}
			}
		}
	}   

	//cout << "No. of upper triangular pruned: " << num_upper_triangular << endl;
	m_vi_Vertices.push_back(m_vi_Edges.size());

	for(int i=0;i < row; i++) 
	{
		m_vi_Edges.insert(m_vi_Edges.end(),nodeList[i].begin(),nodeList[i].end());
		m_vi_Vertices.push_back(m_vi_Edges.size());
	}

	if(b_getValue) 
	{
		for(int i=0;i<row; i++) 
		{
			m_vd_Values.insert(m_vd_Values.end(),valueList[i].begin(),valueList[i].end());
		}
	}

	nodeList.clear();
	valueList.clear();
	// CalculateVertexDegrees();
	return true;

}


}