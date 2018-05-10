/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include "map.h"

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) 
{
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 300;

	default_random_engine gen;
	double std_x, std_y, std_theta;

	std_x = std[0];
	std_y = std[1];
	std_theta = std[2];

	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	for (int i = 0; i < num_particles; ++i) 
	{
		double sample_x, sample_y, sample_theta;

		sample_x = dist_x(gen);
		sample_y = dist_y(gen);
		sample_theta = dist_theta(gen);

		Particle temp;
		temp.x = sample_x;
		temp.y = sample_y;
		temp.theta = sample_theta;
		temp.weight = 1.0;

		particles.push_back(temp);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) 
{
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;

	normal_distribution<double> dist_x(0.0, std_pos[0]);
	normal_distribution<double> dist_y(0.0, std_pos[1]);
	normal_distribution<double> dist_theta(0.0, std_pos[2]);


	vector<Particle> new_particles;
	// cout << yaw_rate << "  YAW RATE" << endl;

	for (int i = 0; i < num_particles; ++i) 
	{
		Particle temp = particles[i];
		
		double current_x = temp.x;
		double current_y = temp.y;
		double current_theta = temp.theta;

		// Got to add here Noises Component
		if (fabs(yaw_rate)>0.001)
		{
			temp.x = current_x + ((velocity / yaw_rate) * (sin(current_theta + (yaw_rate * delta_t))  - sin(current_theta))) + dist_x(gen);
			temp.y = current_y + ((velocity / yaw_rate) * (-cos(current_theta + (yaw_rate * delta_t))  + cos(current_theta))) + dist_y(gen);
			temp.theta = current_theta + (delta_t * yaw_rate) + dist_theta(gen);
			temp.weight = 1.0;
		}
		else
		{
			temp.x = current_x + velocity* cos(current_theta) *delta_t + dist_x(gen);
			temp.y = current_y + velocity* sin(current_theta) *delta_t + dist_y(gen);
			temp.theta = current_theta + (delta_t * yaw_rate) + dist_theta(gen);
			temp.weight = 1.0;
		}

		new_particles.push_back(temp);
	}

	particles = new_particles;

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) 
{
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	for (int i = 0; i < observations.size(); ++i) 
	{
		int current_index = -1000;
		double least_distance = 100000000000;
		for (int j = 0; j < predicted.size(); ++j) 
		{
			double dist = sqrt(pow(predicted[j].x - observations[i].x, 2) + pow(predicted[j].y - observations[i].y, 2));
			if (dist < least_distance)
			{
				current_index = predicted[j].id;
				least_distance = dist;
				// cout << i << " " << j << "  " << dist << endl;
			}
		}
		observations[i].id = current_index; 
	}
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) 
{
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	double denom = 0.0;
	for (int i = 0; i < num_particles; ++i) 
	{
		Particle current_particle = particles[i];
		vector<LandmarkObs> predicted;
		vector<Map::single_landmark_s> landmark_list = map_landmarks.landmark_list;
 
		vector<LandmarkObs> particle_obs = observations;

		for(int k = 0; k < landmark_list.size(); k++)
		{ 
			LandmarkObs single_landmark;
			
			Map::single_landmark_s sls = landmark_list[k];
			single_landmark.id = sls.id_i;
			
			// single_landmark.x = current_particle.x + (cos(current_particle.theta) * sls.x_f) - (sin(current_particle.theta) * sls.y_f);
			// single_landmark.y = current_particle.y + (sin(current_particle.theta) * sls.x_f) + (cos(current_particle.theta) * sls.y_f);
			single_landmark.x = sls.x_f;
			single_landmark.y = sls.y_f;

			double distance = sqrt(pow(current_particle.x - single_landmark.x, 2) + pow(current_particle.y - single_landmark.y, 2));

			if(distance < sensor_range)
			{
				predicted.push_back(single_landmark);
			}
		}

		for(int k = 0; k < particle_obs.size(); k++)
		{ 			 
			double new_x = current_particle.x + (cos(current_particle.theta) * particle_obs[k].x) - (sin(current_particle.theta) * particle_obs[k].y);
			double new_y = current_particle.y + (sin(current_particle.theta) * particle_obs[k].x) + (cos(current_particle.theta) * particle_obs[k].y);
			
			particle_obs[k].y = new_y;
			particle_obs[k].x = new_x;

		}

		dataAssociation(predicted, particle_obs);

		double product_prob = 1.0;
		double pro;
		for(int k = 0; k < particle_obs.size(); k++)
		{
			LandmarkObs single_obs = particle_obs[k];

			for (int p = 0; p < predicted.size() ; p ++)
			{
				LandmarkObs current_landmark = predicted[p];
				if(current_landmark.id == single_obs.id)
				{
					// cout << "OBS == X" << single_obs.x << " LandMark ==   X" << current_landmark.x << endl;
					// cout << "OBS == Y" << single_obs.y << " LandMark ==   Y" << current_landmark.y << endl;
					pro =  ((1/(2 * M_PI * std_landmark[0] * std_landmark[1])) * exp(-1 * ( (pow(single_obs.x - current_landmark.x, 2)/(2 * pow(std_landmark[0], 2))) + (pow(single_obs.y - current_landmark.y, 2)/ (2 * pow(std_landmark[1], 2))))));
					product_prob *= pro; 
					// cout << pro << endl;
					break;
				}
			}
		}

		particles[i].weight = product_prob;
		
		// cout << "Prob == " << product_prob << "  ::  " << endl;

		denom += product_prob;
	}

	for (int i = 0; i < num_particles; ++i) 
	{
		particles[i].weight = particles[i].weight/denom;
		// cout << i << "  +++   " << particles[i].weight << "  ---- " << denom << endl;
	}

}

void ParticleFilter::resample() 
{
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	vector<Particle> new_particles;

	random_device rd;
	mt19937 gen(rd());

	vector<double> weights;

	for(int i = 0; i < num_particles; i++)
	{
		weights.push_back(particles[i].weight);
	}

	discrete_distribution<> d(weights.begin(), weights.end());


	for(int n=0; n < num_particles; ++n) 
	{
		Particle  temp;
        int index = d(gen);
        temp = particles[index];
        new_particles.push_back(temp);
    }

    particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
