/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


import java.io.FileNotFoundException;

/**
 *
 * @author Uporabnik
 */
public class KServerOffline {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws FileNotFoundException {
        Ks k = new Ks(2,13,1);
        k.startGreedy(true);
    }
    
}
