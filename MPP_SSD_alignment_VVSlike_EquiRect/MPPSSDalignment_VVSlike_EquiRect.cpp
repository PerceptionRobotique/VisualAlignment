/*!
 \file MPPSSDalignment_VVSlike_EquiRect.cpp
 \brief Mixture of Photometric Potentials (MPP) SSD for spherical camera orientation estimation (3 DOFs) in VVS-like scheme, exploiting PeR core, core_extended, io, features, estimation and sensor_pose_estimation modules
 * example command line :
 *  ./MPPSSDalignment_VVSlike_EquiRect 3 0.325 ./2021_MPP_SSD_alignment_EquiRect_media/images_subdiv3/ ./2021_MPP_SSD_alignment_EquiRect_media/depthmaps_subdiv3/ 0 0 1 1 ./2021_MPP_SSD_alignment_EquiRect_media/images_subdiv3/maskFull.png 1 1 1 0
 \param subDiv the number of subdivision levels for the spherical image sampling
 \param lambda_g the Gaussian expansion parameter
 \param imDir the directory containing the Equirectangular images to process
 \param depthDir the directory containing the Equirectangular depthmaps to process
 \param iRef the reference image index (in the lexicographical order)
 \param i0 the first image index of the sequence to process
 \param i360 the last image index
 \param iStep the image sequence looping step
 \param Mask the image file of the mask (white pixels are to be considered whereas black pixels are not)
 \param nbTries the number of tested initial guesses for the optimization (the one leading to the lower MPP-SSD is kept)
 \param estimationType selects which estimation type to consider between 0 pure alignment, 4 alignment sequence (later: 1 incremental alignment, 2 incremental alignment with key images)
 \param stabilization if 1, outputs the rotation compensated dualfisheye image
 \param truncGauss truncated Gaussian domain: 1 yes (+ or - 3 lambda_g at most), 0 no (default)
 \param ficPosesInit the text file of initial poses (one pose line per image to process)
 *
 \author Guillaume CARON
 \version 0.1
 \date January 2022
 */

#include <iostream>
#include <iomanip>

#include <per/prStereoModel.h>

#include <per/prStereoModelXML.h>
#include <per/prRegularlySampledCSImageDepth.h>

#include <per/prPhotometricGMS.h>
#include <per/prFeaturesSet.h>

#include <per/prSSDCmp.h>

#include <per/prPoseSphericalEstim.h>

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>

#include <visp/vpImageTools.h>
#include <visp/vpImage.h>
#include <visp/vpImageIo.h>

#include <visp/vpTime.h>

#include <visp/vpDisplayX.h>

#define INTERPTYPE prInterpType::IMAGEPLANE_BILINEAR

//#define VERBOSE

/*!
 * \fn main()
 * \brief Main function of the MPP SSD based spherical alignement
 *
 * 1. Alignment objects initialization, considering the pose estimation of a spherical camera from the feature set of photometric Gaussian mixture 3D samples compared thanks to the SSD
 * 2. Successive computation of the "desired" festures set for every image of the sequence that are used to register the request spherical image considering zero values angles initialization, the optimal angles of the previous image (the request image changes at every iteration), the optimal angles of the previous image (the resquest image changes only if the MPP-SSD error is greater than a threshold)
 * 3. Save the MMP-SSD at optimal poses, optimal poses, processing times and key images numbers to files
 *
 * \return
 *          0 if the program ended without any issue
 *         -2 if no subdivision level is provided
 *         -3 no lambda_g value
 *         -4 no image files directory path
 *         -9 no depthmap files directory path
 *         -5 no reference image file index
 *         -6 no initial file index
 *         -7 no last file index
 *         -8 no image step
 */
int main(int argc, char **argv)
{
    
    //1. Get the number of subdivision levels
    if(argc < 2)
    {
#ifdef VERBOSE
        std::cout << "no subdivision level" << std::endl;
#endif
        return -2;
    }
    unsigned int subdivLevel = atoi(argv[1]);
    
    //Get the lambda_g value
    if(argc < 3)
    {
#ifdef VERBOSE
        std::cout << "no lambda_g" << std::endl;
#endif
        return -3;
    }
    float lambda_g = atof(argv[2]);
    
    //Loading the reference image with respect to which the cost function will be computed
    vpImage<unsigned char> I_req;
    if(argc < 4)
    {
#ifdef VERBOSE
        std::cout << "no image files directory path given" << std::endl;
#endif
        return -4;
    }
    char *imPath = (char *)argv[3];

		//Loading the depthmap of the reference image with respect to which the cost function will be computed
    vpImage<float> depthmap_req;
    if(argc < 5)
    {
#ifdef VERBOSE
        std::cout << "no depthmap files directory path given" << std::endl;
#endif
        return -9;
    }
    char *depthmapPath = (char *)argv[4];

    //Get filename thanks to boost
    char myFilter[1024], myFilter_depth[1024];
    char ext[] = "png";
		char ext_depthmap[] = "pfm";
    if(argc < 6)
    {
#ifdef VERBOSE
        std::cout << "no reference image file number given" << std::endl;
#endif
        return -5;
    }
    unsigned int iRef = atoi(argv[5]); 
    
    if(argc < 7)
    {
#ifdef VERBOSE
        std::cout << "no initial image file number given" << std::endl;
#endif
        return -6;
    }
    unsigned int i0 = atoi(argv[6]); 
    
    if(argc < 8)
    {
#ifdef VERBOSE
        std::cout << "no image files count given" << std::endl;
#endif
        return -7;
    }
    unsigned int i360 = atoi(argv[7]);
    
    if(argc < 9)
    {
#ifdef VERBOSE
        std::cout << "no image sequence step given" << std::endl;
#endif
        return -8;
    }
    unsigned int iStep = atoi(argv[8]);

		//Load and display I_req
    sprintf(myFilter, "%06d.*\\.%s", iRef, ext);
    
    boost::filesystem::path dir(imPath);
    boost::regex my_filter( myFilter );
    std::string name;
    
    for (boost::filesystem::directory_iterator iter(dir),end; iter!=end; ++iter)
    {
        name = iter->path().filename().string();
        if (boost::regex_match(name, my_filter))
        {
            vpImageIo::read(I_req, iter->path().string());
            break;
        }
    }

    vpDisplayX disp;
    disp.init(I_req, 25, 25, "I_req");
    vpDisplay::display(I_req);
    vpDisplay::flush(I_req);
		//vpDisplay::getClick(I_req);

		//Load and display depthmap_req
    sprintf(myFilter_depth, "%06d.*\\.%s", iRef, ext_depthmap); //filter on extension not working
    
    boost::filesystem::path dir_depthmap(depthmapPath);
    boost::regex my_filter_depthmap( myFilter_depth );
    
    for (boost::filesystem::directory_iterator iter(dir_depthmap),end; iter!=end; ++iter)
    {
        name = iter->path().filename().string();
        if (boost::regex_match(name, my_filter_depthmap))
        {
            vpImageIo::readPFM(depthmap_req, iter->path().string());
            break;
        }
    }
		
		//resize first the depthmap to match the image, if necessary
		if( (depthmap_req.getWidth() != I_req.getWidth()) || (depthmap_req.getHeight() != I_req.getHeight()) )
		{
			vpImage<float> depthmap_req_resize;
			vpImageTools::resize(depthmap_req, depthmap_req_resize, I_req.getWidth(), I_req.getHeight(), vpImageTools::INTERPOLATION_LINEAR);
			depthmap_req = depthmap_req_resize;
		}
		
		vpImage<unsigned char> depthmap_req_disp;
		vpImageConvert::convert(depthmap_req, depthmap_req_disp); 	
    vpDisplayX disp_depthmap;

    disp_depthmap.init(depthmap_req_disp, 250, 25, "depthmap_req");
    vpDisplay::display(depthmap_req_disp);
    vpDisplay::flush(depthmap_req_disp);
		//vpDisplay::getClick(depthmap_req_disp);
    
    //load the mask image
    vpImage<unsigned char> Mask;
    if(argc < 10)
    {
#ifdef VERBOSE
        std::cout << "no mask image given" << std::endl;
#endif
        Mask.resize(I_req.getHeight(), I_req.getWidth(), 255);
    }
    else
    {
        try
        {
            vpImageIo::read(Mask, argv[9]);
        }
        catch(vpException e)
        {
            std::cout << "unable to load mask file" << std::endl;
            Mask.resize(I_req.getHeight(), I_req.getWidth(), 255);
        }
    }
    
    //number of tries to define the initial pose rotation r_0
    unsigned int nbTries = 1;
    if(argc < 11)
    {
#ifdef VERBOSE
        std::cout << "no initial number of tries given" << std::endl;
#endif
        //return -9;
    }
    else
        nbTries = atoi(argv[10]);
    
    //estimation type (0: pure alignment ; 1: odometry ; 2 : odometry with key images)
    unsigned int estimationType = 0;
    if(argc < 12)
    {
#ifdef VERBOSE
        std::cout << "no estimation type given" << std::endl;
#endif
        //return -9;
    }
    else
        estimationType = atoi(argv[11]);
    
    //stabilisation de la sequence
    unsigned int stabilisation = 0;
    if(argc < 13)
    {
#ifdef VERBOSE
        std::cout << "no stabilisation parameter given" << std::endl;
#endif
        //return -9;
    }
    else
        stabilisation = atoi(argv[12]);

		//truncated Gaussians
    unsigned int truncGauss = 0;
    if(argc < 14)
    {
#ifdef VERBOSE
        std::cout << "no Gaussian truncature parameter given" << std::endl;
#endif
        //return -9;
    }
    else
        truncGauss = atoi(argv[13]);
    
    //file with initial poses r_0
    bool ficInit = false;
    std::vector<vpPoseVector> v_pv_init;
    if(argc < 15)
    {
#ifdef VERBOSE
        std::cout << "no initial poses file given" << std::endl;
#endif
        //return -9;
    }
    else
    {
        ficInit = true;

        std::ifstream ficPosesInit(argv[14]);
        vpPoseVector r;
        while(!ficPosesInit.eof())
        {
            ficPosesInit >> r[0] >> r[1] >> r[2] >> r[3] >> r[4] >> r[5];
            v_pv_init.push_back(r);
        }
        ficPosesInit.close();
    }

    
    // 2. VVS objects initialization, considering the pose control of a virtual spherical camera from the feature set of photometric Gaussian mixture "3D" samples compared thanks to the SSD

    /*unsigned */int ehaut, elarg;
    ehaut = I_req.getHeight();
    elarg = I_req.getWidth();
    prEquirectangular ecam(elarg*0.5/M_PI, ehaut*0.5/(M_PI*0.5), elarg*0.5, ehaut*0.5);

//    bool dofs[6] = {false, false, false, true, true, true}; //"gyro"
    bool dofs[6] = {true, true, true, true, true, true}; //"6DoF"
//    bool dofs[6] = {false, false, true, false, false, false}; //"vertical" translation DoF only


    //prepare the request spherical image (here, the reference image is always considered as the request in order to compute the 3D transformation that allows to transform it to the target image)
		prRegularlySampledCSImageDepth<unsigned char> IS_req(subdivLevel); //the regularly sample spherical image-depth to be set from the acquired/loaded dual fisheye image
    IS_req.setInterpType(prInterpType::IMAGEPLANE_BILINEAR);
    
    IS_req.buildFromEquiRect(I_req, depthmap_req, ecam, &Mask); 
    IS_req.toAbsZN(); //prepare spherical pixels intensities for the MPP cost function expression constraints
		prRegularlySampledCSImageDepth<float> GS(subdivLevel); //contains every pr3DCartesianPointVec XS_g and does GS_sample.buildFrom(IS_req, XS_g);
    
    prFeaturesSet<prCartesian3DPointVec, prPhotometricGMS<prCartesian3DPointVec>,prRegularlySampledCSImageDepth > fSet_req;
    prPhotometricGMS<prCartesian3DPointVec> GS_sample_req(lambda_g, truncGauss==1);

    bool poseJacobianCompute = true;
    // TODO : calculer en parallele un fSet_req avec lambda_g /= 10 pour les dernières itérations --> précision accrue, sans perdre de temps
    fSet_req.buildFrom(IS_req, GS, GS_sample_req, poseJacobianCompute);
    
    prPhotometricGMS<prCartesian3DPointVec> GS_sample(lambda_g, truncGauss==1);
    std::cout << "nb initial features : " << fSet_req.set.size() << std::endl;
    
    vpDisplayX disp2;
    
    //to save iterations
    std::ostringstream s;
    std::string filename;
    s.str("");
    s.setf(std::ios::right, std::ios::adjustfield);
    s << imPath << "/iter_" << iRef << "_" << i0 << "_" << i360 << ".txt";
    filename = s.str();
      
    //For each image of the dataset
    int nbPass = 0;
    bool clickOut = false;
    unsigned int imNum = i0;
    std::vector<double> err;
    std::vector<vpPoseVector> pv;
    double temps;
    std::vector<double> v_temps;
    std::vector<unsigned int> v_keyImageNum;
    v_temps.reserve((i360-i0)/iStep);
    
    vpPoseVector r, r_to_save;
    vpHomogeneousMatrix key_dMc, dMd_prec;
    
    //de/activate the M-Estimator
    bool robust = false;//true;//
    vpImage<unsigned char> I_des;
    
    //3. Successive computation of the "desired" (target) festures set for every image of the sequence that are used to register the request spherical image considering zero values angles initialization, the optimal angles of the previous image (the request image changes at every iteration), the optimal angles of the previous image (the resquest image changes only if the MPP-SSD error is greater than a threshold)
    prFeaturesSet<prCartesian3DPointVec, prPhotometricGMS<prCartesian3DPointVec>, prRegularlySampledCSImageDepth > fSet_des; //prRegularlySampledCSImageDepth though this depth is not used nor existing
    double seuilErr = 0.0325; //0.015; //0.0077;// // OK for 0,325 only and subdiv3
    while(!clickOut && (imNum <= i360))
    {
        temps = vpTime::measureTimeMs();
        std::cout << "num request image : " << nbPass << std::endl;
        
        switch(estimationType)
        {
            case 0: //pure alignment
            default:
            {
                r.set(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                break;
            }
            case 1: //odometry
            {
                if(nbPass > 0)
                {
                    key_dMc.buildFrom(r_to_save);
                    fSet_req = fSet_des;
                    //virtual_servo.buildFrom(fSet_req); //done below on fSet_des now
                    r.set(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                }
                break;
            }
            case 2: //odometry with key images
            {
                if( (nbPass > 0) && (err[nbPass-1] > seuilErr) )
                {
                    key_dMc.buildFrom(r_to_save);
                    fSet_req = fSet_des; //check si ce n'est pas encore la precedente !
                    //virtual_servo.buildFrom(fSet_req); //done below on fSet_des now
                    v_keyImageNum.push_back(nbPass-1);
                    r.set(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                }
                break;
            }
            case 4: //alignment "sequence" (no reinit)
            {
                break;
            }
        }
        
        //load and display the target image        
        sprintf(myFilter, "%06d.*\\.%s", imNum, ext);
        
        my_filter.set_expression(myFilter);
        
        for (boost::filesystem::directory_iterator iter(dir),end; iter!=end; ++iter)
        {
            name = iter->path().leaf().string();
            if (boost::regex_match(name, my_filter))
            {
                std::cout << iter->path().string() << " loaded" << std::endl;
                vpImageIo::read(I_des, iter->path().string());
                break;
            }
        }
        if(nbPass == 0)
            disp2.init(I_des, 500, 50, "I_des");

        vpDisplay::display(I_des);
        vpDisplay::flush(I_des);

        //I_r - I_des
        s.str("");
        s.setf(std::ios::right, std::ios::adjustfield);
        s << imPath << "/aligned/" << std::setfill('0') << std::setw(6) << imNum << "_diff_0.png";
        filename = s.str();
        vpImage<unsigned char> I_diff;
        vpImageTools::imageDifference(I_req, I_des, I_diff);
        vpImageIo::write(I_diff, filename);
        
        // Desired features set setting from the current image
				prRegularlySampledCSImageDepth<unsigned char> IS_des(subdivLevel);
        IS_des.setInterpType(prInterpType::IMAGEPLANE_BILINEAR);
        IS_des.buildFromEquiRect(I_des, ecam, &Mask); //no depth to mimic the future application case of alignment on a RGB captured image
        IS_des.toAbsZN();
        
        //calculer en parallele un fSet_des avec lambda_g /= 10 pour les dernières itérations --> précision accrue, sans perdre de temps
        fSet_des.buildFrom(IS_des, GS, GS_sample, false); // Goulot !
        std::cout << "nb target features : " << fSet_des.set.size() << std::endl;
        
        // if there is a file provided as initial poses, they are used instead of other strategies
        if(ficInit)
        {
            r = v_pv_init[nbPass];
            std::cout << "r init : " << r.t() << std::endl;
        }
        else
        {
            // trying to select the best initial 3D orientation guess
            if(nbTries > 1)
            {
                vpPoseVector r_best_init;
                double err_min_init = 1e20;
                double angle[3]={0,0,0}, pasAngulaire;
                vpHomogeneousMatrix dMc;
                double err0;
                pasAngulaire = 2.0*M_PI / nbTries;

                unsigned int nbTriesPerDOF[3] = {1,1,1};
                if(dofs[3])
                {
                    nbTriesPerDOF[0] = nbTries;
                    angle[0] = -pasAngulaire*floor(nbTries*0.5);
                }
                if(dofs[4])
                {
                    nbTriesPerDOF[1] = 1;//nbTries;
                    angle[1] = -pasAngulaire*floor(nbTries*0.5);
                }
                if(dofs[5])
                {
                    nbTriesPerDOF[2] = 1;//nbTries;
                    angle[2] = -pasAngulaire*floor(nbTries*0.5);
                }
                for(int iTry0 = 0 ; iTry0 < nbTriesPerDOF[0] ; iTry0++, angle[0]+=pasAngulaire)
                {
                    r[3] = angle[0];
                    if(dofs[4])
                        angle[1] = -pasAngulaire*floor(nbTries*0.5);
                    else
                        angle[1] = 0;
                    for(int iTry1 = 0 ; iTry1 < nbTriesPerDOF[1] ; iTry1++, angle[1]+=pasAngulaire)
                    {
                        r[4] = angle[1];
                        if(dofs[5])
                            angle[2] = -pasAngulaire*floor(nbTries*0.5);
                        else
                            angle[2] = 0;
                        for(int iTry2 = 0 ; iTry2 < nbTriesPerDOF[2] ; iTry2++, angle[2]+=pasAngulaire)
                        {
                            r[5] = angle[2];
                            
                            dMc.buildFrom(r);
                            fSet_req.update(dMc);
                            
                            prSSDCmp<prCartesian3DPointVec, prPhotometricGMS<prCartesian3DPointVec> > errorComputer(fSet_req, fSet_des, robust);
                            prPhotometricGMS<prCartesian3DPointVec> GS_error = errorComputer.getRobustCost();
                            err0 = GS_error.getGMS();
                            
                            if(err0 < err_min_init)
                            {
                                err_min_init = err0;
                                r_best_init = r;
                            }
                        }
                    }
                }
                r = r_best_init;
            }
        }
        
        // register the request feature set over the desired one and save the optimal MPP-SSD
        
        //initialization of (V)VS
    		prPoseSphericalEstim<prFeaturesSet<prCartesian3DPointVec, prPhotometricGMS<prCartesian3DPointVec>, prRegularlySampledCSImageDepth >, 
                         prSSDCmp<prCartesian3DPointVec, prPhotometricGMS<prCartesian3DPointVec> > > virtual_servo;                     

    		virtual_servo.setdof(dofs[0], dofs[1], dofs[2], dofs[3], dofs[4], dofs[5]);   
        
        virtual_servo.buildFrom(fSet_des);
    
    		virtual_servo.initControl(0.7f); // gain, cst_Z
        
				//Control loop
				//3. Successive computation of the current features set for every acquired image that are used to control the camera toward the desired image
		    bool updateSampler = true;
		    bool poseJacobianCompute = true;

				// ----------------------------------------------------------
				unsigned int nbDOF = 6, numDOF, indDOF;
				int iter   = 1;
				vpColVector v6(6), v;
				vpHomogeneousMatrix dMc0, cip1Mci;
				cip1Mci.eye();
				double residual;

				double tms, duree;
				unsigned int nbIter = 0;

				//to display current images at each iteration of the VVS loop
				vpImage<unsigned char> I_r(I_des.getHeight(), I_des.getWidth());
				vpDisplayX disp_r;
    		disp_r.init(I_r, 1000, 25, "I_r");
    		
				do
				{
					tms = vpTime::measureTimeMs();
					std::cout << "--------------------------------------------" << iter++ << std::endl ;

					//later update the spherical image-depth with new pose
					/*
					//update current (named req) features set
					IS_req.buildFromEquiRect(I_req, depthmap_req, ecam, &Mask, &cip1Mci);  //replace with CS image (req?) //IP_cur and I_cur in VisualServoing...
					IS_req.toAbsZN(); 
					fSet_req.updateMeasurement(IS_req, GS, GS_sample_req, poseJacobianCompute, updateSampler); //fSet_cur and IP_cur in VisualServoing...
					*/
					//now, waiting for confirming a buildFrom can safely be done several times on a image-depth:
					prRegularlySampledCSImageDepth<unsigned char> IS_cur(subdivLevel); //the regularly sample spherical image-depth to be set from the acquired/loaded dual fisheye image
					IS_cur.setInterpType(prInterpType::IMAGEPLANE_BILINEAR);
					
					IS_cur.buildFromEquiRect(I_req, depthmap_req, ecam, &Mask, &dMc0); 
					IS_cur.toAbsZN();
		
					//fSet_req.updateMeasurement(IS_cur, GS, GS_sample_req, poseJacobianCompute, updateSampler); //fSet_cur and IP_cur in VisualServoing...
					//now, waiting for updateMeasurement to be validated:
					prFeaturesSet<prCartesian3DPointVec, prPhotometricGMS<prCartesian3DPointVec>,prRegularlySampledCSImageDepth > fSet_cur;
					prPhotometricGMS<prCartesian3DPointVec> GS_sample_cur(lambda_g, truncGauss==1);
					fSet_cur.buildFrom(IS_cur, GS, GS_sample_cur, poseJacobianCompute);

					//Compute control / pose update vector
					residual = 0.5*virtual_servo.control(fSet_cur, v, robust); 
					
					std::cout << "error : " << residual << std::endl;

					//update the DOFs
					indDOF = 0;
					for (numDOF = 0 ; numDOF < nbDOF ; numDOF++)
							if (dofs[numDOF])
							{
							    v6[numDOF] = -v[indDOF];
							    indDOF++;
							}
							else
							    v6[numDOF] = 0;

					std::cout << "v6 : " << v6.t() << std::endl;
					
					//apply computed velocity to dMc0 (transformation matrix from the req/initially rendered image to the target/desired image) with the exp map of se3
					//pose update vector -> transformation matrix
					cip1Mci = vpExponentialMap::direct(v6).inverse();
					dMc0 = cip1Mci*dMc0;
					
					duree = vpTime::measureTimeMs() - tms;
					std::cout << "duration : " << duree <<std::endl;

					//display the current image of the VVS loop
					r.buildFrom(dMc0);	
					dMd_prec.buildFrom(r);
		      r_to_save.buildFrom(dMd_prec*key_dMc);
					I_r = 0;

		      if(stabilisation)
		      {
		          vpPoseVector ir;
		          //r_to_save.set(0,0,0,0,0,0);
		          ir.buildFrom(vpHomogeneousMatrix(r_to_save).inverse());
		          //IS_des.toTwinOmni(I_r, ir, stereoCam, &Mask);
		          IS_des.toEquiRect(I_r, ir, ecam, &Mask);
		      }
		      else
		      {
		      		vpPoseVector ir;
		      		//r.set(0,0,0,0,0,0);
		          //IS_req.toTwinOmni(I_r, r, stereoCam, &Mask);
		          ir.buildFrom(vpHomogeneousMatrix(r_to_save).inverse());
		          IS_req.toEquiRect(I_r, ir, ecam, &Mask);
		      }

					vpDisplay::display(I_r);
		  		vpDisplay::flush(I_r);				
				}
				while((++nbIter < 50) && !vpDisplay::getClick(I_req, false)); //change the VVS loop exit criterion
				
				//save the optimal MPP-SSD
				err.push_back(residual);
				
        v_temps.push_back(vpTime::measureTimeMs()-temps);
        std::cout << "Pass " << nbPass << " time : " << v_temps[nbPass] << " ms" << std::endl;
        
        
        pv.push_back(r_to_save);
        
        std::cout << "Pose optim : " << r.t() << " cum : " << r_to_save.t() << std::endl;
        
        std::cout << "weighted FPP-SSD : " << err[nbPass] << std::endl;
        
        clickOut=vpDisplay::getClick(I_req,false);
        
        
        s.str("");
        s.setf(std::ios::right, std::ios::adjustfield);
        s << imPath << "/aligned/" << std::setfill('0') << std::setw(6) << imNum << ".png";
        filename = s.str();
        vpImageIo::write(I_r, filename);

        //I_r - I_des
        s.str("");
        s.setf(std::ios::right, std::ios::adjustfield);
        s << imPath << "/aligned/" << std::setfill('0') << std::setw(6) << imNum << "_diff_1.png";
        filename = s.str();
        //vpImage<unsigned char> I_diff;
        vpImageTools::imageDifference(I_r, I_des, I_diff);
        vpImageIo::write(I_diff, filename);

        
        
        imNum+=iStep;
        nbPass++;
    }
    
    // 4. Save the MMP-SSD at optimal poses, optimal poses, processing times and key images numbers to files
    
    //save err list to file
    s.str("");
    s.setf(std::ios::right, std::ios::adjustfield);
    s << imPath << "/errors_" << iRef << "_" << i0 << "_" << i360 << ".txt";
    filename = s.str();
    std::ofstream ficerrMin(filename.c_str());
    std::vector<double>::iterator it_err = err.begin();
    for(;it_err != err.end() ; it_err++)
    {
        ficerrMin << *it_err << std::endl;
    }
    ficerrMin.close();
    
    //save poses list to file
    s.str("");
    s.setf(std::ios::right, std::ios::adjustfield);
    s << imPath << "/poses_" << iRef << "_" << i0 << "_" << i360 << ".txt";
    filename = s.str();
    std::ofstream ficPoses(filename.c_str());
    std::vector<vpPoseVector>::iterator it_pv = pv.begin();
    for(;it_pv != pv.end() ; it_pv++)
    {
        ficPoses << it_pv->t() << std::endl;
    }
    ficPoses.close();
    
    //save times list to file
    s.str("");
    s.setf(std::ios::right, std::ios::adjustfield);
    s << imPath << "/time_" << iRef << "_" << i0 << "_" << i360 << ".txt";
    filename = s.str();
    std::ofstream ficTime(filename.c_str());
    std::vector<double>::iterator it_time = v_temps.begin();
    for(;it_time != v_temps.end() ; it_time++)
    {
        ficTime << *it_time << std::endl;
    }
    ficTime.close();
    
    //save key images numbers to file
    s.str("");
    s.setf(std::ios::right, std::ios::adjustfield);
    s << imPath << "/keys_" << iRef << "_" << i0 << "_" << i360 << ".txt";
    filename = s.str();
    std::ofstream ficKeys(filename.c_str());
    std::vector<unsigned int>::iterator it_keys = v_keyImageNum.begin();
    for(;it_keys != v_keyImageNum.end() ; it_keys++)
    {
        ficKeys << *it_keys << std::endl;
    }
    ficKeys.close();
    
	return 0;
}
