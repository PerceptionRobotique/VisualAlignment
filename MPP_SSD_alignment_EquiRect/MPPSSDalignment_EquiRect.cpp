/*!
 \file MPPSSDgyroEstim_EquiRect.cpp
 \brief Mixture of Photometric Potentials (MPP) SSD for spherical camera orientation estimation (3 DOFs), exploiting PeR core, core_extended, io, features, estimation and sensor_pose_estimation modules
 * example command line :
 *  ./MPPSSDalignment_EquiRect 3 0.325 ./2021_MPP_SSD_alignment_EquiRect_media/images_subdiv3/ ./2021_MPP_SSD_alignment_EquiRect_media/depthmaps_subdiv3/ 0 0 1 1 ./2021_MPP_SSD_alignment_EquiRect_media/images_subdiv3/maskFull.png 1 1 1 0
 
 ./MPPSSDalignment_EquiRect 5 0.025 /home/guillaume/Data/FullScan/pcrgbalign/images_subdiv5/ /home/guillaume/Data/FullScan/pcrgbalign/depthmaps_full/ 0 603 603 1 /home/guillaume/Data/FullScan/pcrgbalign/images_subdiv5/mask.png 1 0 0 1
 
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
 \param estimationType selects which estimation type to consider between 0 pure gyro, 1 incremental gyro, 2 incremental gyro with key images
 \param stabilization if 1, outputs the rotation compensated dualfisheye image
 \param truncGauss truncated Gaussian domain: 1 yes (+ or - 3 lambda_g at most), 0 no (default)
 \param ficPosesInit the text file of initial poses (one pose line per image to process)
 *
 \author Guillaume CARON
 \version 0.1
 \date june 2021
 */

#include <iostream>
#include <iomanip>

#include <per/prStereoModel.h>

#include <per/prStereoModelXML.h>
//#include <per/prRegularlySampledCSImage.h>
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
 * \brief Main function of the MPP SSD based spherical orientation estimation
 *
 * 1. Loading a divergent stereovision system made of two fisheye cameras considering the Barreto's model from an XML file got from the MV calibration software
 * 2. Gyro objects initialization, considering the pose estimation of a spherical camera from the feature set of photometric Gaussian mixture 3D samples compared thanks to the SSD
 * 3. Successive computation of the "desired" festures set for every image of the sequence that are used to register the request spherical image considering zero values angles initialization, the optimal angles of the previous image (the request image changes at every iteration), the optimal angles of the previous image (the resquest image changes only if the MPP-SSD error is greater than a threshold)
 * 4. Save the MMP-SSD at optimal poses, optimal poses, processing times and key images numbers to files
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
    float lambda_g = atof(argv[2]);//0.35;//0.035;//
    
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
    unsigned int iRef = atoi(argv[5]); //72
    
    if(argc < 7)
    {
#ifdef VERBOSE
        std::cout << "no initial image file number given" << std::endl;
#endif
        return -6;
    }
    unsigned int i0 = atoi(argv[6]);//1;//0;
    
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

std::cout << I_req.getWidth() << " " << I_req.getHeight() << std::endl;    

		//Load and display I_req
    sprintf(myFilter, "%06d.*\\.%s", iRef, ext);
    
    boost::filesystem::path dir(imPath);
    boost::regex my_filter( myFilter );
    std::string name;
    
    for (boost::filesystem::directory_iterator iter(dir),end; iter!=end; ++iter)
    {
        name = iter->path().filename().string();
std::cout << iter->path().string() << std::endl;
        if (boost::regex_match(name, my_filter))
        {
            vpImageIo::read(I_req, iter->path().string());
            break;
        }
    }

std::cout << I_req.getWidth() << " " << I_req.getHeight() << std::endl;

    vpDisplayX disp;
    disp.init(I_req, 25, 25, "I_req");
    vpDisplay::display(I_req);
    vpDisplay::flush(I_req);
		vpDisplay::getClick(I_req);

std::cout << depthmap_req.getWidth() << " " << depthmap_req.getHeight() << std::endl;

		//Load and display depthmap_req
    sprintf(myFilter_depth, "%06d.*\\.%s", iRef, ext_depthmap); //filter on extension not working
    
    boost::filesystem::path dir_depthmap(depthmapPath);
    boost::regex my_filter_depthmap( myFilter_depth );
    
    for (boost::filesystem::directory_iterator iter(dir_depthmap),end; iter!=end; ++iter)
    {
        name = iter->path().filename().string();
std::cout << iter->path().string() << std::endl;
        if (boost::regex_match(name, my_filter_depthmap))
        {
            vpImageIo::readPFM(depthmap_req, iter->path().string());
            break;
        }
    }

		std::cout << depthmap_req.getWidth() << " " << depthmap_req.getHeight() << std::endl;
		
		if( (depthmap_req.getWidth() != I_req.getWidth()) || (depthmap_req.getHeight() != I_req.getHeight()) )
		{
			vpImage<float> depthmap_req_resize;
			vpImageTools::resize(depthmap_req, depthmap_req_resize, I_req.getWidth(), I_req.getHeight(), vpImageTools::INTERPOLATION_LINEAR);//vpImageTools::INTERPOLATION_AREA);
			depthmap_req = depthmap_req_resize;
		}
		
		std::cout << depthmap_req.getWidth() << " " << depthmap_req.getHeight() << std::endl;
		
		vpImage<unsigned char> depthmap_req_disp;
		vpImageConvert::convert(depthmap_req, depthmap_req_disp); 	
    vpDisplayX disp_depthmap;

    disp_depthmap.init(depthmap_req_disp, 250, 25, "depthmap_req");
    vpDisplay::display(depthmap_req_disp);
    vpDisplay::flush(depthmap_req_disp);

		vpDisplay::getClick(depthmap_req_disp);
    
    //lecture de l'image "masque"
    //Chargement du masque
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
    
    //nombre de coups d'essai pour definir r_0
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
    
    //type d'estimation (0 : gyro pur ; 1 : odometrie ; 2 : odometrie a images cles)
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
    
    //fichier avec les poses initiales r_0
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

    
    // 2. Gyro objects initialization, considering the pose estimation of a spherical camera from the feature set of photometric Gaussian mixture 3D samples compared thanks to the SSD

    /*unsigned */int ehaut, elarg;
    ehaut = I_req.getHeight();
    elarg = I_req.getWidth();
    //nbPixelse = ehaut*elarg;
    prEquirectangular ecam(elarg*0.5/M_PI, ehaut*0.5/(M_PI*0.5), elarg*0.5, ehaut*0.5);
    
    /*
    //En image plane
    prPhotometricGMS<pr2DCartesianPointVec> G_sample;
    pr2DCartesianPointVec u_g;
    G_sample.buildFrom(I_req, u_g, lambda_g);
    */
    
    //En image spherique
    //initialisation de l'estimation d'orientation
    prPoseSphericalEstim<prFeaturesSet<prCartesian3DPointVec, prPhotometricGMS<prCartesian3DPointVec>, prRegularlySampledCSImageDepth >, prSSDCmp<prCartesian3DPointVec, prPhotometricGMS<prCartesian3DPointVec> > > gyro(1e-6);
//    bool dofs[6] = {false, false, false, true, true, true}; //"gyro"
    bool dofs[6] = {true, true, true, true, true, true}; //"6DoF"
//    bool dofs[6] = {true, true, true, false, false, false}; //"3DoF"
//    bool dofs[6] = {false, false, true, false, false, false}; //"1DoF"

    gyro.setdof(dofs[0], dofs[1], dofs[2], dofs[3], dofs[4], dofs[5]);

    //prepare the request spherical image (here, the reference image is always considered as the request in order to compute the rotations that allow to rotate it to the current image)
    //prRegularlySampledCSImageDepth<unsigned char, float> IS_req(subdivLevel); //the regularly sample spherical image-depth to be set from the acquired/loaded dual fisheye image
		prRegularlySampledCSImageDepth<unsigned char> IS_req(subdivLevel); //the regularly sample spherical image-depth to be set from the acquired/loaded dual fisheye image
    IS_req.setInterpType(prInterpType::IMAGEPLANE_BILINEAR);
    
    //IS_req.buildFromTwinOmni(I_req, stereoCam, &Mask); // Goulot !
    IS_req.buildFromEquiRect(I_req, depthmap_req, ecam, &Mask); 
    IS_req.toAbsZN(); //prepare spherical pixels intensities for the MPP cost function expression constraints
		//below: prRegularlySampledCSImage or prRegularlySampledCSImageDepth?
    //prRegularlySampledCSImageDepth<float, float> GS(subdivLevel); //contient tous les pr3DCartesianPointVec XS_g et fera GS_sample.buildFrom(IS_req, XS_g);
		prRegularlySampledCSImageDepth<float> GS(subdivLevel); //contient tous les pr3DCartesianPointVec XS_g et fera GS_sample.buildFrom(IS_req, XS_g);
    
    prFeaturesSet<prCartesian3DPointVec, prPhotometricGMS<prCartesian3DPointVec>,prRegularlySampledCSImageDepth > fSet_req;
    prPhotometricGMS<prCartesian3DPointVec> GS_sample_req(lambda_g, truncGauss==1);

    bool poseJacobianCompute = true;
    // TODO : calculer en parallele un fSet_req avec lambda_g /= 10 pour les dernières itérations --> précision accrue, sans perdre de temps
    fSet_req.buildFrom(IS_req, GS, GS_sample_req, poseJacobianCompute);

    gyro.buildFrom(fSet_req);
    
    prPhotometricGMS<prCartesian3DPointVec> GS_sample(lambda_g, truncGauss==1);
    std::cout << "nb features : " << fSet_req.set.size() << std::endl;
    
    vpDisplayX disp2;
    
    //to save iterations
    std::ostringstream s;
    std::string filename;
    s.str("");
    s.setf(std::ios::right, std::ios::adjustfield);
    s << imPath << "/iter_" << iRef << "_" << i0 << "_" << i360 << ".txt";
    filename = s.str();
    gyro.startSaveIterations((char *)filename.c_str());
    
    //Pour chaque image du dataset
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
    
    //activate the M-Estimator
    bool robust = false;//true;//
    vpImage<unsigned char> I_des;
    vpImage<unsigned char> I_diff;
    
    //3. Successive computation of the "desired" festures set for every image of the sequence that are used to register the request spherical image considering zero values angles initialization, the optimal angles of the previous image (the request image changes at every iteration), the optimal angles of the previous image (the resquest image changes only if the MPP-SSD error is greater than a threshold)
    //double angle = -177.5*M_PI/180.;
    prFeaturesSet<prCartesian3DPointVec, prPhotometricGMS<prCartesian3DPointVec>, prRegularlySampledCSImageDepth > fSet_des; //prRegularlySampledCSImageDepth though this depth is not used nor existing
    double seuilErr = 0.0325; //0.015; //0.0077;// // OK pour 0,325 seul et subdiv3
    while(!clickOut && (imNum <= i360))
    {
        temps = vpTime::measureTimeMs();
        std::cout << "num request image : " << nbPass << std::endl;
        
        switch(estimationType)
        {
            case 0: //gyro pur
            default:
            {
                r.set(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                break;
            }
            case 1: //odometrie
            {
                if(nbPass > 0)
                {
                    key_dMc.buildFrom(r_to_save);
                    fSet_req = fSet_des;
                    gyro.buildFrom(fSet_req);
                    r.set(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                }
                break;
            }
            case 2: //odometrie a images cles
            {
                if( (nbPass > 0) && (err[nbPass-1] > seuilErr) )
                {
                    key_dMc.buildFrom(r_to_save);
                    fSet_req = fSet_des; //check si ce n'est pas encore la precedente !
                    gyro.buildFrom(fSet_req);
                    v_keyImageNum.push_back(nbPass-1);
                    r.set(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
                }
                break;
            }
            case 4: //gyro "sequence" (pas de reinit)
            {
                break;
            }
        }
        
        sprintf(myFilter, "%06d.*\\.%s", imNum, ext);
        
        my_filter.set_expression(myFilter);
        
        for (boost::filesystem::directory_iterator iter(dir),end; iter!=end; ++iter)
        {
            name = iter->path().leaf().string();
            if (boost::regex_match(name, my_filter))
            {
                std::cout << iter->path().string() << " loaded" << std::endl;
                vpImageIo::read(I_des, iter->path().string());
                
								vpImageTools::imageDifference(I_req, I_des, I_diff);
								s.str("");
								s.setf(std::ios::right, std::ios::adjustfield);
								s << imPath << "/aligned/" << std::setfill('0') << std::setw(6) << imNum << "I_diff_init.png";
								filename = s.str();
								vpImageIo::write(I_diff, filename);
								
                break;
            }
        }
        if(nbPass == 0)
            disp2.init(I_des, 500, 50, "I_des");

        vpDisplay::display(I_des);
        vpDisplay::flush(I_des);
        
        // Desired feature set setting from the current image
        //prRegularlySampledCSImageDepth<unsigned char, float> IS_des(subdivLevel);
				prRegularlySampledCSImageDepth<unsigned char> IS_des(subdivLevel);
        IS_des.setInterpType(prInterpType::IMAGEPLANE_BILINEAR);
        //IS_des.buildFromTwinOmni(I_des, stereoCam, &Mask);
        IS_des.buildFromEquiRect(I_des, ecam, &Mask); //no depth
        IS_des.toAbsZN();
        
        //calculer en parallele un fSet_des avec lambda_g /= 10 pour les dernières itérations --> précision accrue, sans perdre de temps
        fSet_des.buildFrom(IS_des, GS, GS_sample, false); // Goulot !
        std::cout << "nb features : " << fSet_des.set.size() << std::endl;
        
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
        err.push_back(gyro.track(fSet_des, r, robust));
    
        v_temps.push_back(vpTime::measureTimeMs()-temps);
        std::cout << "Pass " << nbPass << " time : " << v_temps[nbPass] << " ms" << std::endl;
        
        dMd_prec.buildFrom(r);
        r_to_save.buildFrom(dMd_prec*key_dMc);
        pv.push_back(r_to_save);
        
        std::cout << "Pose optim : " << r.t() << " cum : " << r_to_save.t() << std::endl;
        
        std::cout << "weighted FPP-SSD : " << err[nbPass] << std::endl;
        
        clickOut=vpDisplay::getClick(I_req,false);
        
        vpImage<unsigned char> I_r(I_des.getHeight(), I_des.getWidth());
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
        		//r.set(0,0,-1,0,0,0);
            //IS_req.toTwinOmni(I_r, r, stereoCam, &Mask);
            IS_req.toEquiRect(I_r, r, ecam, &Mask);
        }
        s.str("");
        s.setf(std::ios::right, std::ios::adjustfield);
        s << imPath << "/aligned/" << std::setfill('0') << std::setw(6) << imNum << ".png";
        filename = s.str();
        vpImageIo::write(I_r, filename);
        
        vpImageTools::imageDifference(I_r, I_des, I_diff);
        s.str("");
        s.setf(std::ios::right, std::ios::adjustfield);
        s << imPath << "/aligned/" << std::setfill('0') << std::setw(6) << imNum << "I_diff_final.png";
        filename = s.str();
        vpImageIo::write(I_diff, filename);
        
        imNum+=iStep;
        nbPass++;
        //angle += 2.5*M_PI/180.;
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
