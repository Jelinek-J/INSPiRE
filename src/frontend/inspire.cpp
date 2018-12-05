// inspire.cpp : Defines the entry point for the console application.
// Last error id: 18

#include "../backend/index.h"
#include "../backend/features.h"
#include "../backend/subgraphs.h"
#include "../backend/fingerprints.h"
#include "../backend/mine.h"
#include "../backend/classify.h"
#include "../backend/predict.h"
#include "../common/filesystem.h"
#include "../common/string.h"
#include <iostream>
#include <filesystem>

//#define TESTING
//#define MAKE

static const std::string CONFIGURATION_FILE = "config";
static const std::string INDEX_FILE = "residues.ind";
static const std::string COORDINATES_FILE = "coordinates.tur";
static const std::string INTERFACES_FILE = "interfaces.tur";
static const std::string FEATURES_FILE = "features.tur";
static const std::string RADIUSES_FILE = "radiuses.rus";
static const std::string AMINOACID_FILE = "aminoacid.moc";
static const std::string NODES_FILE = "nodes.sup";
static const std::string EDGES_FILE = "edges.sup";
static const std::string QUERY_FILE = "fingerprints.fit";
static const std::string FINGERPRINTS_FILE = "settings.json";
static const std::string MINED_FILE = "mined.med";
static const std::string STATISTICS_FILE = "ratios.sas";
static const std::string PREDICTION_FILE = "prediction.pec";
static const double DEFAULT_DISTANCE = 0.5;

// Initialize dictionary file for aminoacids
static void create_aminoacid_file(std::string path) {
  std::ofstream aminoacids(path);
  aminoacids << "00A\t5\n02K\tA\n02L\tN\n02O\tA\n02Y\tA\n033\tV\n037\tP\n03Y\tC\n04U\tP\n04V\tP\n05N\tP\n07O\tC\n0A0\tD\n0A1\tY\n0A2\tK\n0A8\tC\n0A9\tF\n0AA\tV\n0AB\tV\n0AC\tG\n";
  aminoacids << "0AD\t2\n0AF\tW\n0AG\tL\n0AH\tS\n0AK\tD\n0AM\t0\n0AP\t1\n0AR\tR\n0AU\t4\n0AV\t0\n0BN\tF\n0CS\tA\n0E5\tT\n0EA\tY\n0FL\tA\n0LF\tP\n0NC\tA\n0QL\tC\n0R8\t1\n0SP\t0\n";
  aminoacids << "0TD\tD\n0UH\t2\n0UO\tW\n0WZ\tY\n0X9\tR\n0Y8\tP\n102\t7\n10C\t6\n11Q\tP\n11W\tE\n125\t9\n126\t9\n127\t9\n12A\t5\n12L\tP\n12X\tP\n12Y\tP\n143\tC\n18M\t7\n18Q\t4\n";
  aminoacids << "1AC\tA\n1AP\t0\n1CC\t1\n1FC\t1\n1L1\tA\n1MA\t5\n1MG\t7\n1OP\tY\n1PA\tF\n1PI\tA\n1RN\t9\n1SC\t6\n1TQ\tW\n1TY\tY\n1X6\tS\n200\tF\n23F\tF\n23G\t7\n23P\tA\n26A\t5\n";
  aminoacids << "26B\tT\n28X\tT\n2AG\tA\n2AR\t0\n2AT\t3\n2AU\t9\n2BT\t3\n2BU\t0\n2CO\tC\n2DA\t0\n2DT\t3\n2EG\t7\n2FM\tM\n2GT\t3\n2GX\tF\n2HF\tH\n2JG\tS\n2JV\t2\n2KK\tK\n2KP\tK\n";
  aminoacids << "2LT\tY\n2LU\tL\n2MA\t5\n2MG\t7\n2ML\tL\n2MR\tR\n2MT\tP\n2MU\t9\n2NT\t3\n2OM\t9\n2OR\tR\n2OT\t3\n2P0\tP\n2PR\t2\n2QZ\tT\n2R3\tY\n2RA\tA\n2RX\tS\n2SG\t7\n2SO\tH\n";
  aminoacids << "2ST\t3\n2TY\tY\n2VA\tV\n2XA\tC\n2ZC\tS\n30F\tU\n30V\tC\n31H\t5\n31M\t5\n31Q\tC\n33S\tF\n33W\tA\n34E\tV\n3AH\tH\n3AU\t9\n3BY\tP\n3CF\tF\n3CT\tY\n3DA\t0\n3GA\tA\n";
  aminoacids << "3GL\tE\n3MD\tD\n3ME\t9\n3MU\t9\n3MY\tY\n3NF\tY\n3O3\tE\n3PX\tP\n3QN\tK\n3TD\t9\n3TT\tP\n3WS\tA\n3WX\tP\n3X9\tC\n3XH\tG\n3YM\tY\n3ZH\tH\n41H\tF\n41Q\tN\n42Y\tS\n";
  aminoacids << "432\tS\n45F\tP\n47C\t1\n4AF\tF\n4AK\tK\n4AR\tR\n4AW\tW\n4BF\tY\n4CF\tF\n4CY\tM\n4D4\tR\n4DP\tW\n4FB\tP\n4FW\tW\n4GJ\tC\n4HH\tS\n4HJ\tS\n4HL\tY\n4HT\tW\n4IN\tW\n";
  aminoacids << "4J4\tC\n4KY\tP\n4L0\tP\n4LZ\tY\n4MM\tM\n4N7\tP\n4N8\tP\n4N9\tP\n4OC\t6\n4OU\tF\n4OV\tS\n4OZ\tS\n4PC\t1\n4PD\t1\n4PE\t1\n4PH\tF\n4PQ\tW\n4SC\t1\n4SJ\tF\n4SU\t9\n";
  aminoacids << "4U3\t1\n4U7\tA\n4WQ\tA\n51T\tY\n54C\tW\n56A\tH\n5AA\t0\n5AB\tA\n5AT\t3\n5BU\t9\n5CG\t2\n5CM\t1\n5CR\tF\n5CS\tC\n5CT\tK\n5CW\tW\n5FA\t5\n5FC\t1\n5FQ\tA\n5FU\t9\n";
  aminoacids << "5GG\tK\n5GM\tI\n5HC\t1\n5HM\t6\n5HP\tE\n5HT\t3\n5HU\t4\n5IC\t6\n5IT\t3\n5IU\t4\n5JP\tS\n5MC\t6\n5MU\t9\n5MW\tK\n5NC\t1\n5OH\tA\n5OW\tK\n5PC\t1\n5PG\tG\n5PY\t3\n";
  aminoacids << "5SE\t4\n5VV\tN\n5XU\tA\n60F\tC\n63G\t2\n63H\t2\n64T\t3\n66D\tI\n68Z\t2\n6BR\tT\n6CL\tK\n6CT\t3\n6CV\tA\n6CW\tW\n6DN\tK\n6FK\t2\n6G4\tK\n6GL\tA\n6HA\t0\n6HB\t0\n";
  aminoacids << "6HC\t1\n6HG\t2\n6HN\tK\n6HT\t3\n6IA\t5\n6M6\tC\n6MA\t5\n6MC\t5\n6MP\t5\n6MT\t5\n6MZ\t5\n6NW\t5\n6OG\t2\n6OO\t6\n6PO\t2\n6V1\tC\n6WK\tC\n6Y9\tP\n70U\t9\n73C\tS\n";
  aminoacids << "73N\tR\n73O\tY\n73P\tK\n73W\t6\n74P\tK\n75B\t9\n77Y\t4\n7AT\t5\n7BG\t2\n7DA\t0\n7GU\t2\n7JA\tI\n7MG\t7\n7N8\tF\n823\tN\n85F\tC\n85Y\t9\n8AA\t7\n8AG\t2\n8AN\t5\n";
  aminoacids << "8AY\tA\n8BA\t0\n8DT\t4\n8FG\t2\n8JB\tC\n8LJ\tP\n8MG\t2\n8OG\t2\n8OS\t7\n8PY\t2\n8RE\tK\n8RJ\t8\n8RO\t1\n8SP\tS\n8WY\tL\n999\tD\n9DN\tN\n9DS\tG\n9E7\tK\n9NE\tE\n";
  aminoacids << "9NF\tF\n9NR\tR\n9NV\tV\n9QV\t9\n9WV\tA\nA\t5\nA23\t5\nA2L\t5\nA2M\t5\nA34\t0\nA35\t0\nA38\t0\nA39\t5\nA3A\t0\nA3P\t5\nA40\t0\nA43\t0\nA44\t5\nA47\t0\nA5L\t0\n";
  aminoacids << "A5M\t6\nA5N\tN\nA5O\t5\nA6A\t5\nA6C\t6\nA6G\t7\nA6U\t9\nA7E\t5\nA8E\tV\nA9D\tS\nA9Z\t5\nAA3\tA\nAA4\tA\nAAR\tR\nABA\tA\nABR\t0\nABS\t0\nACL\tR\nAD2\t0\nADI\t5\n";
  aminoacids << "ADP\t5\nAEI\tD\nAET\t5\nAF2\t0\nAFA\tN\nAFG\t2\nAGM\tR\nAGQ\tY\nAGT\tC\nAHB\tN\nAHL\tR\nAHO\tA\nAHP\tA\nAIB\tA\nAKL\tD\nAKZ\tD\nALA\tA\nALC\tA\nALM\tA\nALN\tA\n";
  aminoacids << "ALO\tT\nALS\tA\nALT\tA\nALV\tA\nALY\tK\nAMD\t5\nAMO\t5\nAN6\tL\nAN8\tA\nAP7\t5\nAPI\tK\nAPK\tK\nAR2\tR\nAR4\tE\nAR7\tR\nARG\tR\nARM\tR\nARO\tR\nAS\t0\nAS7\tN\n";
  aminoacids << "ASA\tD\nASB\tD\nASI\tD\nASK\tD\nASL\tD\nASN\tN\nASP\tD\nASQ\tD\nASX\tB\nATD\t3\nATL\t3\nATM\t3\nAVC\t5\nAVJ\tH\nAYA\tA\nAZH\tA\nAZK\tK\nAZS\tS\nAZY\tY\nB1F\tF\n";
  aminoacids << "B27\tT\nB2A\tA\nB2F\tF\nB2I\tI\nB2V\tV\nB3A\tA\nB3D\tD\nB3E\tE\nB3K\tK\nB3U\tH\nB3X\tN\nB3Y\tY\nB7C\t1\nB8H\t9\nB8K\t7\nB8Q\t6\nB8T\t6\nB8W\t7\nB9B\t7\nB9H\t6\n";
  aminoacids << "BB6\tC\nBB7\tC\nBB8\tF\nBB9\tC\nBBC\tC\nBCS\tC\nBCX\tC\nBFD\tD\nBG1\tS\nBGH\t7\nBGM\t2\nBH2\tD\nBHD\tD\nBIF\tF\nBIU\tI\nBL2\tL\nBLE\tL\nBLY\tK\nBMT\tT\nBNN\tF\n";
  aminoacids << "BOE\t3\nBOR\tR\nBP5\tA\nBPE\tC\nBRU\t4\nBSE\tS\nBTA\tL\nBTC\tC\nBTK\tK\nBTR\tW\nBUC\tC\nBUG\tV\nBVP\t4\nBWV\tR\nBYR\tY\nC\t6\nC1J\tR\nC1S\tC\nC1T\tC\nC1X\tK\n";
  aminoacids << "C22\tA\nC25\t6\nC2L\t6\nC2S\t1\nC31\t6\nC32\t1\nC34\t1\nC36\t1\nC37\t1\nC38\t1\nC3Y\tC\nC42\t1\nC43\t6\nC45\t1\nC46\t1\nC49\t1\nC4G\tR\nC4R\tC\nC4S\t1\nC5C\tC\n";
  aminoacids << "C5L\t6\nC67\tR\nC6C\tC\nC6D\tR\nC6G\t2\nC7R\t1\nC7S\t1\nCAF\tC\nCAR\t1\nCAS\tC\nCAY\tC\nCB2\t1\nCBR\t1\nCBV\t6\nCCC\t6\nCCS\tC\nCDW\t1\nCE7\tN\nCEA\tC\nCFL\t1\n";
  aminoacids << "CFZ\t1\nCG1\t7\nCG6\tC\nCGA\tE\nCGU\tE\nCGV\tC\nCH\t6\nCHP\tG\nCIR\tR\nCLE\tL\nCLG\tK\nCLH\tK\nCME\tC\nCMH\tC\nCML\tC\nCMR\t1\nCMT\tC\nCNU\t9\nCP1\t1\nCR5\tG\n";
  aminoacids << "CS0\tC\nCS1\tC\nCS3\tC\nCS4\tC\nCSA\tC\nCSB\tC\nCSD\tC\nCSE\tC\nCSF\t6\nCSJ\tC\nCSL\t1\nCSO\tC\nCSP\tC\nCSR\tC\nCSS\tC\nCSU\tC\nCSW\tC\nCSX\tC\nCSZ\tC\nCTE\tW\n";
  aminoacids << "CTG\t3\nCTH\tT\nCWD\tA\nCWR\tS\nCX2\t1\nCXM\tM\nCY0\tC\nCY1\tC\nCY3\tC\nCY4\tC\nCYA\tC\nCYD\tC\nCYF\tC\nCYG\tC\nCYJ\tK\nCYM\tC\nCYQ\tC\nCYR\tC\nCYS\tC\nCYW\tC\n";
  aminoacids << "CZ2\tC\nCZS\tA\nCZZ\tC\nD00\t1\nD2T\tD\nD3T\t3\nD4M\t3\nDA\t0\nDA2\tR\nDAB\tA\nDAH\tF\nDBS\tS\nDBU\tT\nDBY\tY\nDBZ\tA\nDC\t1\nDC2\tC\nDCG\t2\nDCT\t1\nDDE\tH\n";
  aminoacids << "DDG\t2\nDDN\t4\nDDZ\tA\nDFC\t1\nDFG\t2\nDG\t2\nDG8\t2\nDGI\t2\nDGP\t2\nDHA\tS\nDHN\tV\nDHU\t9\nDI7\tY\nDIR\tR\nDLS\tK\nDM0\tK\nDMH\tN\nDMK\tD\nDNL\tK\nDNP\tA\n";
  aminoacids << "DNR\t1\nDNS\tK\nDNW\tA\nDOC\t1\nDOH\tD\nDON\tL\nDP1\tR\nDPB\t3\nDPL\tP\nDPP\tA\nDPQ\tY\nDRM\t4\nDRT\t3\nDT\t3\nDU\t4\nDUZ\t4\nDYA\tD\nDYJ\tP\nDYS\tC\nDZM\t0\n";
  aminoacids << "E\t0\nE0Y\tP\nE1X\t0\nE3C\t6\nE6G\t7\nE7G\t7\nE9V\tH\nECC\tQ\nECX\tC\nEDA\t0\nEFC\tC\nEFG\t2\nEHG\t2\nEHP\tF\nEIT\t3\nELY\tK\nEME\tE\nEPM\tM\nEPQ\tQ\nESB\tY\n";
  aminoacids << "ESC\tM\nEUP\tT\nEXA\tK\nEXC\t1\nEXY\tL\nF2F\tF\nF2T\t9\nF3H\t3\nF3N\t5\nF4H\t3\nFA2\t0\nFAK\tK\nFB5\tA\nFB6\tA\nFC0\tF\nFCL\tF\nFDG\t2\nFDL\tK\nFFM\tC\nFGL\tG\n";
  aminoacids << "FGP\tS\nFH7\tK\nFHL\tK\nFHO\tK\nFHU\t9\nFIO\tR\nFLA\tA\nFLE\tL\nFLT\tY\nFME\tM\nFMG\t2\nFNU\t9\nFOE\tC\nFOX\t2\nFP9\tP\nFPK\tP\nFT6\tW\nFTR\tW\nFTY\tY\nFVA\tV\n";
  aminoacids << "FY2\tY\nFY3\tY\nFZN\tK\nG\t7\nG01\tE\nG1G\t7\nG25\t7\nG2L\t7\nG2S\t2\nG31\t2\nG32\t2\nG33\t2\nG36\t2\nG38\t2\nG42\t2\nG46\t7\nG47\t2\nG48\t7\nG49\t2\nG7M\t7\n";
  aminoacids << "G8M\tE\nGAO\t7\nGAU\tE\nGCK\t1\nGDO\t7\nGDP\t7\nGDR\t2\nGEE\tG\nGF2\t2\nGFL\t2\nGFT\tS\nGH3\t7\nGHC\tE\nGHG\tQ\nGHW\tE\nGL3\tG\nGLH\tQ\nGLJ\tE\nGLK\tE\nGLN\tQ\n";
  aminoacids << "GLQ\tE\nGLU\tE\nGLX\tZ\nGLY\tG\nGLZ\tG\nGMA\tE\nGME\tE\nGMS\t2\nGMU\t4\nGN7\t2\nGNC\tQ\nGNG\t7\nGOM\t7\nGPL\tK\nGRB\t7\nGS\t2\nGSC\tG\nGSR\t2\nGSS\t2\nGSU\tE\n";
  aminoacids << "GT9\tC\nGTP\t7\nGVL\tS\nGX1\t2\nH14\tF\nH1D\tM\nH2U\t9\nH5M\tP\nHAC\tA\nHAR\tR\nHBN\tH\nHCM\tC\nHDP\t4\nHEU\t4\nHGY\tG\nHHI\tH\nHIA\tH\nHIC\tH\nHIP\tH\nHIQ\tH\n";
  aminoacids << "HIS\tH\nHIX\tA\nHL2\tL\nHL5\tL\nHLU\tL\nHLY\tK\nHMR\tR\nHN0\t2\nHN1\t2\nHNC\tC\nHOX\tF\nHPC\tF\nHPE\tF\nHPH\tF\nHPQ\tF\nHQA\tA\nHR7\tR\nHRG\tR\nHRP\tW\nHS8\tH\n";
  aminoacids << "HS9\tH\nHSE\tS\nHSK\tH\nHSL\tS\nHSO\tH\nHSV\tH\nHT7\tW\nHTI\tC\nHTR\tW\nHV5\tA\nHVA\tV\nHY3\tP\nHYI\tM\nHYP\tP\nHZP\tP\nI2M\tI\nI4G\tG\nI4U\t9\nI58\tK\nI5C\t1\n";
  aminoacids << "IAM\tA\nIAR\tR\nIC\t6\nICY\tC\nIEL\tK\nIG\t7\nIGL\tG\nIGU\t2\nIIL\tI\nILE\tI\nILG\tE\nILM\tI\nILX\tI\nILY\tK\nIMC\t1\nIML\tI\nIMP\t7\nIOR\tR\nIPG\tG\nIT1\tK\n";
  aminoacids << "IU\t9\nIYR\tY\nIZO\tM\nJDT\t3\nJJJ\tC\nJJK\tC\nJJL\tC\nK1R\tC\nK5L\tS\nKAG\t7\nKBE\tK\nKCR\tK\nKCX\tK\nKEO\tK\nKFP\tK\nKGC\tK\nKNB\tA\nKOR\tM\nKPF\tK\nKPI\tK\n";
  aminoacids << "KPY\tK\nKST\tK\nKYN\tW\nKYQ\tK\nL3O\tL\nL5P\tK\nLA2\tK\nLAA\tD\nLAL\tA\nLAY\tL\nLBY\tK\nLC\t6\nLCA\t5\nLCG\t2\nLCK\tK\nLCX\tK\nLDG\t2\nLDH\tK\nLE1\tV\nLED\tL\n";
  aminoacids << "LEF\tL\nLEH\tL\nLEM\tL\nLEN\tL\nLET\tK\nLEU\tL\nLEX\tL\nLG\t7\nLGP\t2\nLGY\tK\nLHU\t9\nLLO\tK\nLLP\tK\nLLY\tK\nLLZ\tK\nLME\tE\nLMF\tK\nLMQ\tQ\nLNE\tL\nLNM\tL\n";
  aminoacids << "LP6\tK\nLPD\tP\nLPG\tG\nLPS\tS\nLRK\tK\nLSH\t3\nLSO\tK\nLST\t3\nLTR\tW\nLVG\tG\nLVN\tV\nLWY\tP\nLYF\tK\nLYK\tK\nLYM\tK\nLYN\tK\nLYP\tK\nLYR\tK\nLYS\tK\nLYU\tK\n";
  aminoacids << "LYX\tK\nLYZ\tK\nM0H\tC\nM1G\t2\nM2G\t7\nM2L\tK\nM2S\tM\nM30\tG\nM3L\tK\nM3R\tK\nM4C\t6\nM5M\t6\nM7A\t5\nMA\tA\nMA6\t5\nMA7\t0\nMAA\tA\nMAD\t5\nMAI\tR\nMBQ\tY\n";
  aminoacids << "MC1\tS\nMCL\tK\nMCS\tC\nMCY\t1\nMD3\tC\nMD5\tC\nMD6\tG\nMDF\tY\nME0\tM\nME6\t1\nMEA\tF\nMEG\tE\nMEN\tN\nMEP\t9\nMEQ\tQ\nMET\tM\nMEU\tG\nMFN\tE\nMFO\t2\nMG1\t2\n";
  aminoacids << "MGG\tR\nMGN\tQ\nMGQ\t5\nMGT\t7\nMGV\t7\nMGY\tG\nMH1\tH\nMH6\tS\nMHG\t7\nMHL\tL\nMHO\tM\nMHS\tH\nMHU\tF\nMIA\t5\nMIR\tS\nMIS\tS\nMK8\tL\nML3\tK\nMLE\tL\nMLL\tL\n";
  aminoacids << "MLY\tK\nMLZ\tK\nMME\tM\nMMO\tR\nMMT\t3\nMNL\tL\nMNU\t9\nMNV\tV\nMP8\tP\nMPQ\tG\nMRG\t2\nMSA\tG\nMSE\tM\nMSL\tM\nMSO\tM\nMT2\tM\nMTR\t3\nMTU\t5\nMTY\tY\nMVA\tV\n";
  aminoacids << "MYK\tK\nMYN\tR\nN10\tS\nN5M\t6\nN6G\t7\nN79\t5\nN7P\tP\nN80\tP\nNA8\tA\nNAL\tA\nNAM\tA\nNBQ\tY\nNC1\tS\nNCB\tA\nNCU\t1\nNDN\t4\nNDU\t4\nNEM\tH\nNEP\tH\nNFA\tF\n";
  aminoacids << "NIY\tY\nNLB\tL\nNLE\tL\nNLN\tL\nNLO\tL\nNLP\tL\nNLQ\tQ\nNLW\tL\nNLY\tG\nNMC\tG\nNMM\tR\nNMS\t3\nNMT\t3\nNNH\tR\nNOT\tL\nNPH\tC\nNPI\tA\nNTR\tY\nNTT\t3\nNTY\tY\n";
  aminoacids << "NVA\tV\nNWD\tA\nNYB\tC\nNYS\tC\nNZC\tT\nNZH\tH\nO2G\t7\nOAR\tR\nOAS\tS\nOBS\tK\nOCS\tC\nOCY\tC\nOGX\t2\nOHI\tH\nOHS\tD\nOHU\t4\nOLD\tH\nOLT\tT\nOLZ\tS\nOMC\t6\n";
  aminoacids << "OMG\t7\nOMH\tS\nOMT\tM\nOMU\t9\nOMX\tY\nOMY\tY\nONE\t9\nONH\tA\nORN\tA\nORQ\tR\nOSE\tS\nOTH\tT\nOXX\tD\nOYL\tH\nP\t2\nP1L\tC\nP2Q\tY\nP2T\t3\nP2U\t4\nP2Y\tP\n";
  aminoacids << "P3Q\tY\nP4U\t9\nP5P\t5\nP7G\t7\nP9S\tC\nPAQ\tY\nPAS\tD\nPAT\tW\nPBB\tC\nPBF\tF\nPCA\tE\nPCC\tP\nPCS\tF\nPDU\t4\nPE1\tK\nPEC\tC\nPF5\tF\nPFF\tF\nPG1\tS\nPG7\t2\n";
  aminoacids << "PGN\t2\nPGP\t7\nPGY\tG\nPH6\tP\nPHA\tF\nPHD\tD\nPHE\tF\nPHI\tF\nPHL\tF\nPHM\tF\nPKR\tP\nPLJ\tP\nPM3\tF\nPMT\t6\nPOM\tP\nPPN\tF\nPPU\t5\nPPW\t2\nPR3\tC\nPR4\tP\n";
  aminoacids << "PR5\t5\nPR7\tP\nPR9\tP\nPRJ\tP\nPRK\tK\nPRN\t0\nPRO\tP\nPRS\tP\nPRV\tG\nPSA\tF\nPSH\tH\nPST\t3\nPSU\t9\nPSW\tU\nPTH\tY\nPTM\tY\nPTR\tY\nPU\t5\nPVH\tH\nPVX\t1\n";
  aminoacids << "PXU\tP\nPYA\tA\nPYH\tK\nPYL\tO\nPYO\t9\nPYX\tC\nPZG\t2\nQBT\t8\nQCS\tC\nQIL\tI\nQMM\tQ\nQPA\tC\nQPH\tF\nQUO\t7\nR\t0\nR1A\tC\nR4K\tW\nRBD\t0\nRDG\t2\nRE0\tW\n";
  aminoacids << "RE3\tW\nRGL\tR\nRIA\t5\nRMP\t0\nRPC\t6\nRPI\tR\nRSQ\t6\nRT\t8\nRT0\tP\nRUS\t9\nRVX\tS\nRZ4\tS\nS12\tS\nS1H\tS\nS2C\tC\nS2M\t3\nS2P\tA\nS4A\t0\nS4C\t6\nS4G\t2\n";
  aminoacids << "S4U\t9\nS6G\t2\nSAC\tS\nSAH\tC\nSAR\tG\nSBG\tS\nSBL\tS\nSC\t1\nSCH\tC\nSCS\tC\nSCY\tC\nSD4\tN\nSDB\tS\nSDE\t0\nSDG\t2\nSDH\t2\nSDP\tS\nSE7\tU\nSEB\tS\nSEC\tU\n";
  aminoacids << "SEE\tS\nSEG\tA\nSEL\tS\nSEM\tS\nSEN\tS\nSEP\tS\nSER\tS\nSET\tS\nSGB\tS\nSHC\tC\nSHP\tG\nSHR\tK\nSIB\tC\nSKH\tK\nSLL\tK\nSLZ\tK\nSMC\tC\nSME\tM\nSMF\tF\nSMP\t0\n";
  aminoacids << "SMT\t3\nSNC\tC\nSNN\tN\nSOC\tC\nSOY\tS\nSPT\t3\nSRA\t5\nSRZ\tS\nSSU\t9\nSTY\tY\nSUN\tS\nSUR\t9\nSVA\tS\nSVV\tS\nSVW\tS\nSVX\tS\nSVY\tS\nSVZ\tS\nSXE\tS\nSYS\tC\n";
  aminoacids << "T\t8\nT0I\tY\nT11\tF\nT23\t8\nT2S\t8\nT31\t9\nT32\t3\nT36\t3\nT37\t3\nT38\t3\nT39\t3\nT3P\t3\nT41\t3\nT48\t3\nT49\t3\nT4S\t3\nT5O\t4\nT5S\t3\nT64\t3\nT6A\t5\n";
  aminoacids << "T8L\tT\nTA3\t3\nTAF\t3\nTAV\tD\nTBG\tV\nTBM\tT\nTBN\t5\nTC1\t1\nTCP\t3\nTCQ\tY\nTCR\tW\nTCY\t0\nTDY\t3\nTED\t3\nTEF\tF\nTFE\t3\nTFF\t3\nTFO\t0\nTFQ\tF\nTFT\t3\n";
  aminoacids << "TGP\t2\nTH5\tT\nTH6\tT\nTHC\tT\nTHR\tT\nTHZ\tR\nTIH\tA\nTIS\tS\nTLC\t3\nTLN\t4\nTLY\tK\nTMB\tT\nTMD\tT\nTNB\tC\nTNR\tS\nTNY\tT\nTOQ\tW\nTOX\tW\nTP1\t3\nTPC\t1\n";
  aminoacids << "TPG\t7\nTPJ\tP\nTPK\tP\nTPL\tW\nTPO\tT\nTPQ\tY\nTQI\tW\nTQQ\tW\nTQZ\tC\nTRF\tW\nTRG\tK\nTRN\tW\nTRO\tW\nTRP\tW\nTRQ\tW\nTRW\tW\nTRX\tW\nTRY\tW\nTS9\tI\nTSP\t8\n";
  aminoacids << "TSY\tC\nTTD\t3\nTTI\t4\nTTM\t3\nTTQ\tW\nTTS\tY\nTXD\t5\nTXP\t5\nTXY\tY\nTY1\tY\nTY2\tY\nTY3\tY\nTY5\tY\nTY8\tY\nTY9\tY\nTYB\tY\nTYI\tY\nTYJ\tY\nTYN\tY\nTYO\tY\n";
  aminoacids << "TYQ\tY\nTYR\tY\nTYS\tY\nTYT\tY\nTYW\tY\nTYY\tY\nU\t9\nU25\t9\nU2L\t9\nU2N\t4\nU2P\t9\nU2X\tY\nU31\t9\nU33\t4\nU34\t9\nU36\t9\nU37\t9\nU3X\tF\nU8U\t9\nUAR\t9\n";
  aminoacids << "UBB\t9\nUBD\t9\nUBI\t4\nUBR\t4\nUCL\t4\nUD5\t9\nUF0\tS\nUF2\t4\nUFR\t4\nUFT\t4\nUGY\tG\nUM1\tA\nUM2\tA\nUMA\tA\nUMO\t9\nUMS\t4\nUMX\t4\nUNK\tX\nUOX\tU\nUPE\t4\n";
  aminoacids << "UPS\t4\nUPV\t9\nUR3\t9\nURD\t9\nURX\t4\nUS1\t4\nUS2\t4\nUS3\t3\nUS5\t9\nUSM\t4\nUVX\t4\nUZR\t9\nV3L\t5\nVAD\tV\nVAF\tV\nVAH\tV\nVAI\tV\nVAL\tV\nVB1\tK\nVH0\tP\n";
  aminoacids << "VR0\tR\nWFP\tF\nWLU\tL\nWPA\tF\nWRP\tW\nWVL\tV\nX\t2\nX2W\tE\nXAD\t0\nXAL\t0\nXCL\t1\nXCN\tC\nXCR\t1\nXCT\t1\nXCY\t1\nXGL\t2\nXGR\t2\nXGU\t2\nXLE\tJ\nXOK\tK\n";
  aminoacids << "XPB\t2\nXPL\tO\nXPR\tP\nXSN\tN\nXTF\t3\nXTH\t3\nXTL\t3\nXTR\t3\nXTS\t7\nXUA\t0\nXUG\t2\nXW1\tA\nXX1\tK\nXYC\tA\nY\t0\nYCM\tC\nYCO\t1\nYG\t7\nYOF\tY\nYPR\tP\n";
  aminoacids << "YPZ\tY\nYTH\tT\nYYG\t7\nZ\t1\nZ01\tA\nZ3E\tT\nZ70\tH\nZAD\t5\nZBC\t6\nZBU\t9\nZBZ\tC\nZCL\tF\nZCY\t6\nZDU\t4\nZGU\t7\nZTH\t8\nZU0\tT\nZYJ\tP\nZYK\tP\nZZD\tC\nZZJ\tA\n";
  aminoacids.flush();
  aminoacids.close();
}

// Prints an information about this program
static void help() {
  std::cout << "Help\n\n";

  std::cout << "Create index of residues in the knowledge-base.\n";
  std::cout << "inspire -h\t\t\tPrints this help.";
  std::cout << "inspire (([-s( <structures>)*]|[-S <structure>]|[-x<temporary_directory>]|[-k<knowledge_base>]) )+\n";
  std::cout << "        [-(b|c|bc|w)] [-i[<radiuses> [<distance>]]] [-F (-a[<aminoacids>]|-t)* -f] ((-n|-e)[c[<int>]|d[<double>]|[e[<double>[-<int>]])*\n";
  std::cout << "        [-g<format>] (([-p<threads>]|[-n<number>]|[-o<center>]) )* [-j<exclude>] [-l<threshold>] [-q<output>]\n";
  std::cout << "inspire (([-s( <structures>)*]|[-S <structure>]|[-x<temporary_directory>]|[-k<knowledge_base>]) )+ -m\n";
  std::cout << "        [-(b|c|bc|w)] [-i[<radiuses> [<distance>]]] [-F (-a[<aminoacids>]|-t)* -f] ((-n|-e)[c[<int>]|d[<double>]|[e[<double>[-<int>]])*\n";
  std::cout << "        [-g<format>]\n";

  std::cout << "\t-s\tList of PDB files to be processed, if <structures> is a directory, it is traversed recursively.\n";
  std::cout << "\t\t<structures> path cannot starts with '-', in that case, -S switcher must be used instead.\n";
  std::cout << "\t-x\tWhere should be stored temporary files. If specified, temporary files are not deleted after processing.\n";
  std::cout << "\t-k\tWhere knowledge-base is (prediction mode)/ should be (construction mode) stored.\n";
  std::cout << "\t-m\tSwitch to construction mode (the default mode is the prediction mode).\n\n";

  std::cout << "\t-cb\tAll biomolecules, models and crystallographic transformations are used\n";
  std::cout << "\t-b\tAll biomolecules and models, but only the first crystallographic transformation are used.\n";
  std::cout << "\t-c\tAll crystallographic transformations, but only the first biomolecule and model are used.\n";
  std::cout << "\t-w\tIgnore both biomolecules and crystallographic transformation, use all chains as they are.\n";
  std::cout << "\t\tIf not switch is typed, only the first biomolecule, model and crystallographic transformations are used.\n\n";

  std::cout << "\t-i\tSpecify a path to a file radiuses of atoms and optionally also the maximal distance between them to be classified as an interface.\n";
  std::cout << "\t-F\tSome features should be extracted and their switchers follows until '-f' occurs.\n";
  std::cout << "\t-a\tAminoacid type should be extracted. Optionally, there can be specified a custom transformation from 3-letters codes to 1-letters codes.\n";
  std::cout << "\t-t\tTemperature of residues should be extracted.\n\n";

  std::cout << "\t-n\tDefine type of neighborhood to be used to generate fingerprints.\n";
  std::cout << "\t-e\tDefine, what residues are connected by edge.\n";
  std::cout << "\tc\tK-nearest residues should be taken.\n";
  std::cout << "\td\tAll residues within <double> Angstroems radius should be taken.\n";
  std::cout << "\te\tAll residues remoted at most <int> 'edges' should be taken, where residues are connected by an edge if they are at most <double> Angstroems distant.\n";
  std::cout << "\t\tDiference between e.g. -nd12.0 and -ne6.0-2 is that in the second case, subgraphs are always connected if edges are defined as -ed6.0.\n\n";

  std::cout << "\t-g\tDefine a custom format of fingerprints.\n";
  std::cout << "\t-p\tChange number of parallel processes to mine similar fingerprints. WARNING: This function is compiller-dependent.\n";
  std::cout << "\t-n\tAt least how many most similar fingerprints should be taken.\n";
  std::cout << "\t-o\tWhat features should be taken into account on central residues. Multiple features must be separatet by directory separator.\n";
  std::cout << "\t-j\tWhat fingerprints should be ignored during mining.\n\n";

  std::cout << "\t-l\tUse linear classifier i.e. if x samples have label 'a' and y samples have label 'b', it will labeled as 'a', if threshold < x/(x+y).\n";
  std::cout << "\t-q\tPrint output in the file <output> instead of on the standart output.\n\n";
}

int main(int argc, const char** argv) {
#ifdef TESTING
  // Errors log for testing reasons
  std::ofstream log;

#ifdef MAKE
  argc = 7;
  const char* arg[] = {argv[0], "-s", "C:\\Inspire\\pdb\\kb\\", "-xC:\\Inspire\\test\\construction\\", "-kC:\\Inspire\\test\\fingerprints", "-m", "-b"};
  log.open("C:\\Inspire\\test\\error-make.log", std::ofstream::app);
#else
  argc = 8;
  const char* arg[] = {argv[0], "-s", "C:\\Inspire\\gvin\\", "-xC:\\Inspire\\basic2\\gvin2\\", "-kC:\\Inspire\\basic2\\fingerprints\\", "-p4", "-jC:\\Inspire\\basic2\\gvin\\related.exc", "-qC:\\Inspire\\basic2\\gvin2\\output"};
  log.open("C:\\Inspire\\basic2\\gvin2\\error-predict.log", std::ofstream::app);
#endif
  argv = arg;
#endif // TESTING

  if (argc > 1 && common::string::starts_with(argv[1], "-h")) {
    help();
    return 0;
  }


  try {
    // Whether the app run in prediction, or construction mode
    bool predict = true;
    // Complexes to precess
    std::vector<std::string> structures;
    // Where the knowledge base is/ should be store(d)
    std::string knowledge_base = "";
    // Filename where read/ write configuration
    std::string config_name;
    // Filestream where read/ write configuration
    std::fstream config_file;
    // Where to store temporary files
    std::string temp_dir;
    // Whether to delete temporary files and directory
    bool delete_tmp = true;
    // Index in argv
    size_t argv_index = 0;

#pragma region Initialization
    {
      while (++argv_index < argc) {
        if (strlen(argv[argv_index]) > 1 && argv[argv_index][0] == '-') {
          if (argv[argv_index][1] == 's') {
            while (++argv_index < argc && (strlen(argv[argv_index]) == 0 || argv[argv_index][0] != '-')) {
              structures.push_back(argv[argv_index]);
            }
            --argv_index;
          } else if (argv[argv_index][1] == 'S') {
            if (++argv_index == argc) {
              std::cerr << "Missing file/ directory specifier after -S switcher";
              return 2224;
            }
            structures.push_back(common::filesystem::enclose_directory_name(argv[argv_index]));
          } else if (argv[argv_index][1] == 'x') {
            temp_dir = common::filesystem::enclose_directory_name(std::string(argv[argv_index]).substr(2));
            delete_tmp = false;
          } else if (argv[argv_index][1] == 'k') {
            knowledge_base = common::filesystem::enclose_directory_name(std::string(argv[argv_index]).substr(2));
          } else if (argv[argv_index][1] == 'm') {
            predict = false;
            ++argv_index;
            break;
          } else {
            break;
          }
        } else {
          std::cerr << "Unexpected switcher '" << argv[argv_index] << "' on position #" << argv_index << ".";
          return 1;
        }
      }

      config_name = knowledge_base + CONFIGURATION_FILE;
      if (predict) {
        if (!common::filesystem::is_directory(knowledge_base)) {
          std::cerr << "Knowledge base '" << knowledge_base << "' must be an existing directory during prediction mode.";
          return 3;
        }
        if (!common::filesystem::is_regular_file(config_name)) {
          std::cerr << "Configuration file '" << config_name << "' does not exists in the input directory or is not a regular file.";
          return 4;
        }
        config_file.open(config_name, std::fstream::in);
      } else {
        if (!(common::filesystem::exists(knowledge_base) || common::filesystem::create_directory_recursive(knowledge_base))) {
          std::cerr << "It is not possible to create knowledge base directory '" << knowledge_base << "'.\n"
            << "  It can occur e.g. if it is an invalid name or some prefix of the path corresponds to an existing file that is not a directory.";
          return 5;
        }
        config_file.open(config_name, std::fstream::out | std::fstream::trunc);
      }
      if (!config_file.is_open()) {
        std::cerr << "It is not possible to open configuration file '" << config_name << "'.";
        return 6;
      }
      if (delete_tmp) {
        for (size_t i = 0; i < 10; i++) {
          std::string tmp = common::filesystem::unique_name();
          if (common::filesystem::create_directory_recursive(tmp)) {
            temp_dir = common::filesystem::enclose_directory_name(tmp);
            break;
          }
        }
        if (temp_dir.size() == 0) {
          std::cerr << "It is not possible to create unique temporary directory.\n  Please try it again or explicitly specify the temporary directory.";
          return 7;
        }
      } else {
        if (!(common::filesystem::is_directory(temp_dir) || common::filesystem::create_directory_recursive(temp_dir))) {
          std::cerr << "It is not possible to create temporary directory '" << temp_dir << "'.\n"
            << "  It can occur e.g. if it is an invalid name or some prefix of the path corresponds to an existing file that is not a directory.";
          return 8;
        }
      }
    }
#pragma endregion Initialization

    std::string index_name = temp_dir + INDEX_FILE;
    inspire::backend::ProteinIterator* it;
    inspire::backend::BasicFilter* filter = new inspire::backend::BasicFilter();
#pragma region Index
    {
      if (predict) {
        std::string bin;
        if (!std::getline(config_file, bin)) {
          std::cerr << "Unexpected end of the configuration file '" << config_name << "'.";
          return 9;
        }
        it = new inspire::backend::ExplicitIterator();
      } else {
        std::string arg_it;
        if (argv_index < argc) {
          arg_it = argv[argv_index];
        }
        if (arg_it == "-c") {
          it = new inspire::backend::FirstModelCrystallographicIterator();
          ++argv_index;
        } else if (arg_it == "-b") {
          it = new inspire::backend::BiomoleculesIterator();
          ++argv_index;
        } else if (arg_it == "-bc") {
          it = new inspire::backend::AllExceptAltLocIterator();
          ++argv_index;
        } else if (arg_it == "-w") {
          it = new inspire::backend::ExplicitIterator();
          ++argv_index;
        } else {
          arg_it = "";
          it = new inspire::backend::FirstModelIterator();
        }
        config_file << arg_it << std::endl;
      }
      inspire::backend::Indexer indexer(index_name, it, filter);
      for (auto structures_it = structures.begin(); structures_it != structures.end(); ++structures_it) {
        indexer.process(*structures_it);
      }
    }
#pragma endregion Index

    std::string coordinates_name = temp_dir + COORDINATES_FILE;
    std::string interfaces_name = knowledge_base + INTERFACES_FILE;
    std::string features_name = temp_dir + FEATURES_FILE;
#pragma region Features
    {
      // Load structure files
      inspire::backend::Features extractor(index_name, it, filter);
      for (auto structures_it = structures.begin(); structures_it != structures.end(); ++structures_it) {
        try {
          extractor.index_pdb(*structures_it);
        } catch (const common::exception::TitledException& e) {
          std::cerr << "ERROR: " << e.what() << std::endl;
        } catch (const std::exception& e) {
          std::cerr << "ERROR: " << e.what() << std::endl;
        } catch (...) {
          std::cerr << "UNKNOWN ERROR" << std::endl;
        }
      }

      // Extract coordinates
      std::vector<inspire::backend::Feature*> features;
      features.push_back(new inspire::backend::CoordinateFeature(it));
      extractor.extract_features(coordinates_name, features);
      delete features.back();

      // Identify interfaces
      if (predict) {
        std::string bin;
        if (!std::getline(config_file, bin)) {
          std::cerr << "Unexpected end of the configuration file '" << config_name << "'.";
          return 18;
        }
      } else {
        features.clear();
        std::string radiuses_name = knowledge_base + RADIUSES_FILE;
        if (predict) {
          std::string arg_dist;
          if (!std::getline(config_file, arg_dist)) {
            std::cerr << "Unexpected end of the configuration file '" << config_name << "'.";
            return 10;
          }
          features.push_back(new inspire::backend::InterfaceLargeFeature(it, radiuses_name, std::stod(arg_dist)));
        } else {
          if (argv_index < argc && common::string::starts_with(argv[argv_index], "-i")) {
            std::string arg(argv[argv_index]);
            size_t space = arg.find_last_of(' ');
            std::string path;
            double dist = DEFAULT_DISTANCE;
            if (space != arg.npos && space != arg.size()-1 && arg[arg.size()-1] != '"' && arg[arg.size()-1] != '\'') {
              path = arg.substr(2, space-2);
              dist = std::stod(arg.substr(space+1));
            } else {
              path = arg.substr(2);
            }
            features.push_back(new inspire::backend::InterfaceLargeFeature(it, path, dist));
            config_file << dist << std::endl;
            common::filesystem::copy(path, radiuses_name);
            ++argv_index;
          } else {
            std::ofstream radiuses(radiuses_name);
            radiuses << "H\t-1.20\nD\t-1.15\nC\t1.70\nN\t1.55\nO\t1.52\nS\t1.80\nP\t1.80\nX\t1.535\n";
            radiuses.flush();
            radiuses.close();
            features.push_back(new inspire::backend::InterfaceLargeFeature(it, radiuses_name, DEFAULT_DISTANCE));
            config_file << DEFAULT_DISTANCE << std::endl;
          }
        }
        extractor.extract_features(interfaces_name, features);
        delete features.back();
      }

      // Extract optional features
      features.clear();
      std::string aminoacid_name = knowledge_base + AMINOACID_FILE;
      if (predict) {
        std::string line;
        while (std::getline(config_file, line) && line != "--") {
          if (line == "-a") {
            features.push_back(new inspire::backend::BasicAminoacidFeature(it, aminoacid_name));
          } else if (line == "-t") {
            features.push_back(new inspire::backend::TemperatureFeature(it));
          } else {
            std::cerr << "Unexpected line '" << line << "' in the configuration file '" << config_name << "'.";
            return 11;
          }
        }
      } else {
        if (argv_index < argc && common::string::starts_with(argv[argv_index], "-F")) {
          while (++argv_index < argc && argv[argv_index] != "-f") {
            if (common::string::starts_with(argv[argv_index], "-a")) {
              if (strlen(argv[argv_index]) > 2) {
                std::string path = std::string(argv[argv_index]).substr(2);
                features.push_back(new inspire::backend::BasicAminoacidFeature(it, path));
                common::filesystem::copy(path, aminoacid_name);
              } else {
                create_aminoacid_file(aminoacid_name);
                features.push_back(new inspire::backend::BasicAminoacidFeature(it, aminoacid_name));
              }
              config_file <<  "-a" << std::endl;
            } else if (argv[argv_index] == std::string("-t")) {
              features.push_back(new inspire::backend::TemperatureFeature(it));
              config_file <<  "-t" << std::endl;
            } else {
              std::cerr << "Unexpected parameter '" << argv[argv_index] << "' within the features section.";
              return 12;
            }
          }
          if (argv_index < argc) {
            ++argv_index;
          }
        } else {
          create_aminoacid_file(aminoacid_name);
          features.push_back(new inspire::backend::BasicAminoacidFeature(it, aminoacid_name));
          config_file <<  "-a" << std::endl;
        }
        config_file << "--" << std::endl;
      }
      extractor.extract_features(features_name, features);
      for (auto features_it = features.begin(); features_it != features.end(); ++features_it) {
        delete *features_it;
      }
    }
#pragma endregion Features
    delete filter;
    delete it;

    std::string nodes_name = temp_dir + NODES_FILE;
    std::string edges_name = temp_dir + EDGES_FILE;
#pragma region Subgraphs
    {
      std::vector<std::string> nodes_edges;
      if (predict) {
        for (size_t i = 0; i < 2; i++) {
          std::string tmp;
          if (!std::getline(config_file, tmp)) {
            std::cerr << "Unexpected end of the configuration file '" << config_name << "'.";
            return 13;
          }
          nodes_edges.push_back(tmp);
        }
      } else {
        nodes_edges.push_back("c12");
        nodes_edges.push_back("d6");
        while (argv_index < argc && strlen(argv[argv_index]) > 1 && argv[argv_index][0] == '-') {
          if (argv[argv_index][1] == 'v') {
            nodes_edges[0] = std::string(argv[argv_index]).substr(2);
          } else if (argv[argv_index][1] == 'e') {
            nodes_edges[1] = std::string(argv[argv_index]).substr(2);
          } else {
            break;
          }
          ++argv_index;
        }
      }

      std::string prefix = temp_dir + "subgraphs-";
      inspire::backend::ChainSubgraphs subgraphs(index_name, coordinates_name, prefix);
      for (size_t i = 0; i < nodes_edges.size(); i++) {
        if (nodes_edges[i].empty()) {
          if (!predict) {
            config_file << "c0" << std::endl;
          }
          nodes_edges[i] = subgraphs.k_nearest(0);
        } else {
          switch ((nodes_edges[i])[0]) {
            case 'c':
              if (nodes_edges[i].size() > 1) {
                if (!predict) {
                  config_file << nodes_edges[i] << std::endl;
                }
                nodes_edges[i] = subgraphs.k_nearest(stoi(std::string(nodes_edges[i]).substr(1)));
              } else {
                if (i == 0) {
                  if (!predict) {
                    config_file << "c15" << std::endl;
                  }
                  nodes_edges[i] = subgraphs.k_nearest(15);
                } else {
                  if (!predict) {
                    config_file << "c6" << std::endl;
                  }
                  nodes_edges[i] = subgraphs.k_nearest(6);
                }
              }
              break;
            case 'd':
              if (nodes_edges[i].size() > 1) {
                if (!predict) {
                  config_file << nodes_edges[i] << std::endl;
                }
                nodes_edges[i] = subgraphs.distance_limit(stof(std::string(nodes_edges[i]).substr(1)));
              } else {
                if (i == 0) {
                  if (!predict) {
                    config_file << "d12" << std::endl;
                  }
                  nodes_edges[i] = subgraphs.distance_limit(12.0);
                } else {
                  if (!predict) {
                    config_file << "d6" << std::endl;
                  }
                  nodes_edges[i] = subgraphs.distance_limit(6.0);
                }
              }
              break;
            case 'e':
              if (nodes_edges[i].size() > 1) {
                std::string parameters(nodes_edges[i]);
                size_t index = parameters.find('-', 2);
                if (index == parameters.npos) {
                  if (i == 0) {
                    if (!predict) {
                      config_file << nodes_edges[i] << "-2" << std::endl;
                    }
                    nodes_edges[i] = subgraphs.edge_limit(stof(parameters.substr(1)), 2);
                  } else {
                    if (!predict) {
                      config_file << nodes_edges[i] << "-1" << std::endl;
                    }
                    nodes_edges[i] = subgraphs.edge_limit(stof(parameters.substr(1)), 1);
                  }
                } else if (index == parameters.size()-1) {
                  if (i == 0) {
                    if (!predict) {
                      config_file << nodes_edges[i] << '2' << std::endl;
                    }
                    nodes_edges[i] = subgraphs.edge_limit(stof(parameters.substr(1, index-1)), 2);
                  } else {
                    if (!predict) {
                      config_file << nodes_edges[i] << '1' << std::endl;
                    }
                    nodes_edges[i] = subgraphs.edge_limit(stof(parameters.substr(1, index-1)), 1);
                  }
                } else {
                  if (!predict) {
                    config_file << nodes_edges[i] << std::endl;
                  }
                  nodes_edges[i] = subgraphs.edge_limit(stof(parameters.substr(1, index-1)), stoi(parameters.substr(index+1)));
                }
              } else {
                if (i == 0) {
                  if (!predict) {
                    config_file << "e6-2" << std::endl;
                  }
                  nodes_edges[i] = subgraphs.edge_limit(6.0, 2);
                } else {
                  if (!predict) {
                    config_file << "e6-1" << std::endl;
                  }
                  nodes_edges[i] = subgraphs.edge_limit(6.0, 1);
                }
              }
              break;
            default:
              std::cerr << "Unknown " << (i == 0 ? "nodes" : "edges") << " identifier '-" << (i == 0 ? 'n' : 'e') << nodes_edges[i] << "'.";
              return 14;
              break;
          }
        }
      }
      subgraphs.extract_subgraphs();
      subgraphs.clear();
      if (nodes_edges[0] == nodes_edges[1]) {
        common::filesystem::copy(nodes_edges[0], nodes_name);
      } else {
        common::filesystem::move(nodes_edges[0], nodes_name);
      }
      common::filesystem::move(nodes_edges[1], edges_name);
    }
#pragma endregion Subgraphs

    std::string query_name = temp_dir + QUERY_FILE;
#pragma region Fingerprints
    {
      inspire::backend::FingerprintWriter fingerprints(index_name, nodes_name);
      fingerprints.add_features(features_name);
      std::string fingerprints_name = knowledge_base + FINGERPRINTS_FILE;
      if (predict) {
        fingerprints.process(fingerprints_name, edges_name, query_name, inspire::backend::FingerprintFormat::Text);
      } else {
        if (argv_index < argc && argv[argv_index] == std::string("-g")) {
          common::filesystem::copy(std::string(argv[argv_index++]).substr(2), fingerprints_name);
        } else {
          std::ofstream settings(fingerprints_name);
          settings << "{\"fingerprint\":{\"size\":1023,\"edge\":[{\"type\":\"distance\",\"size\":4}],\"vertex\": [{\"type\": \"mapping\",\"property\": \"aminoacid\",\"map\": {"
                   << "\"A\": 1,\"B\": 2,\"C\": 3,\"D\": 4,\"E\": 5,\"F\": 6,\"G\": 7,\"H\": 8,\"I\": 9,\"J\": 10,\"K\": 11,\"L\": 12,\"M\": 13,\"N\":14,\"O\": 15,"
                   << "\"P\": 16,\"Q\": 17,\"R\": 18,\"S\": 19,\"T\": 20,\"U\": 21,\"V\": 22,\"W\": 23,\"X\": 0,\"Y\": 24,\"Z\": 25,\"\": 26},\"size\": 5}]}}\n";
          settings.flush();
          settings.close();
        }
        fingerprints.process(fingerprints_name, edges_name, knowledge_base, inspire::backend::FingerprintFormat::Binary);
      }
    }
#pragma endregion Fingerprints

    std::string mined_name = temp_dir + MINED_FILE;
#pragma region Mining
    if (predict) {
      std::set<std::string> filters;
      int threads = 1;
      int limit = 1;
      for (; argv_index < argc && strlen(argv[argv_index]) > 1 && argv[argv_index][0] == '-'; ++argv_index) {
        if (argv[argv_index][1] == 'p') {
          if (strlen(argv[argv_index]) == 2) {
            std::cerr << "Missing integer in number of parallel threads switcher.";
            return 15;
          }
          threads = std::stoi(std::string(argv[argv_index]).substr(2));
        } else if (argv[argv_index][1] == 'n') {
          if (strlen(argv[argv_index]) == 2) {
            std::cerr << "Missing integer in switcher saying how many the most similar fingerprints should be taken.";
            return 16;
          }
          limit = std::stoi(argv[argv_index]);
        } else if (argv[argv_index][1] == 'o') {
          std::stringstream parts(argv[argv_index]);
          char bin;
          parts >> bin >> bin;
          std::string part;
          while (std::getline(parts, part, common::filesystem::directory_separator)) {
            filters.insert(part);
          }
        } else {
          break;
        }
      }
      inspire::backend::Mine mine(knowledge_base, filters, threads, limit);
      for (; argv_index < argc && common::string::starts_with(argv[argv_index], "-j"); ++argv_index) {
        mine.load_excludes(std::string(argv[argv_index]).substr(2));
      }
      mine.select(query_name, mined_name);
    }
#pragma endregion Mining

    std::string statistics_name = temp_dir + STATISTICS_FILE;
#pragma region Classify
    if (predict) {
      inspire::backend::Classifier classifier(interfaces_name);
      classifier.classify(mined_name, statistics_name);
    }
#pragma endregion Classify

#pragma region Predict
    if (predict) {
      inspire::backend::Predictor* predictor = nullptr;
      if (argv_index < argc && common::string::starts_with(argv[argv_index], "-l")) {
        if (strlen(argv[argv_index]) == 2) {
          predictor = new inspire::backend::FractionalPredictor("0.5175");
        } else {
          predictor = new inspire::backend::FractionalPredictor(std::string(argv[argv_index]).substr(2));
        }
        ++argv_index;
      }
      if (predictor == nullptr) {
        predictor = new inspire::backend::FractionalPredictor("0.5175");
      }
      if (argv_index < argc && common::string::starts_with(argv[argv_index], "-q")) {
        std::string prediction_name = std::string(argv[argv_index++]).substr(2);
        if (statistics_name.empty() || statistics_name.back() == common::filesystem::directory_separator) {
          prediction_name += STATISTICS_FILE;
        } else if (common::filesystem::is_directory(statistics_name)) {
          prediction_name += common::filesystem::directory_separator;
          prediction_name += STATISTICS_FILE;
        }
        predictor->predict(statistics_name, prediction_name);
      } else {
        std::string prediction_name = temp_dir + PREDICTION_FILE;
        predictor->predict(statistics_name, prediction_name);
        std::ifstream prediction(prediction_name);
        std::string line;
        while (std::getline(prediction, line)) {
          std::cout << line << '\n';
        }
        prediction.close();
      }
      
      delete predictor;
    }
#pragma endregion Predict

#pragma region Finalization
    if (argv_index < argc) {
      std::cerr << "Some argument" << ((argc - argv_index > 1) ? "s were" : " was") << " not used: '";
      if (argc - argv_index > 1) {
        std::cerr << argv[argv_index];
        for (size_t i = argv_index; ++i < argc-1;) {
          std::cerr << "', '" << argv[i];
        }
        std::cerr << "' and '";
      }
      std::cerr << argv[argc-1] << "'. It can occur if the first of them is not recognized.";
      return 17;
    }
    //TODO: Check whether everything works correctly prior deleting something
    if (false && delete_tmp) {
      common::filesystem::remove_recursively(temp_dir);
    }
#pragma endregion Finalization
  } catch (const common::exception::TitledException& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
#ifdef TESTING
    log << "    ERROR: " << e.what() << std::endl;
#endif // TESTING
    return -1;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
#ifdef TESTING
    log << "    ERROR: " << e.what() << std::endl;
#endif // TESTING
    return -2;
  } catch (...) {
    std::cerr << "UNKNOWN ERROR" << std::endl;
#ifdef TESTING
    log << "    UNKNOWN ERROR" << std::endl;
#endif // TESTING
    return -3;
  }

  return 0;
}
