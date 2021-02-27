import org.jlab.jnp.hipo4.io.HipoReader
import org.jlab.jnp.hipo4.data.Bank
import org.jlab.jnp.hipo4.data.Event
import org.jlab.jnp.hipo4.data.SchemaFactory
import org.jlab.clas.physics.LorentzVector
import uconn.utils.pid.stefan.ElectronCandidate
import uconn.utils.pid.stefan.ElectronCandidate.Cut

def banknames = ['RUN::config', 'REC::Particle','REC::Calorimeter','REC::Cherenkov','REC::Traj']

args.each{fname->
  def reader = new HipoReader()
  reader.open(fname)
  def event = new Event()
  def factory = reader.getSchemaFactory()
  def banks = banknames.collect{new Bank(factory.getSchema(it))}

  while(reader.hasNext()) {
    reader.nextEvent(event)
    banks.each{event.read(it)}

    if(banks.every()) {
      def (runb,partb,calb,ccb,trajb) = banks

      def ipart = 0
      def candidate = ElectronCandidate.getElectronCandidate(ipart, partb, calb, ccb, trajb)

      //test for all cuts specified in ElectronCandidate class
      if(candidate.iselectron()) {
        println("good electron is found")
        def ele = candidate.getLorentzVector()
        println("electron momentum = "+ele.p())
      }

      //specify the list of cuts you want to apply
      if(candidate.iselectron(Cut.ELE_PID, Cut.CC_NPHE, Cut.DC_VERTEX)) {
        println("good electron is found")
        def ele = candidate.getLorentzVector()
        println("electron momentum = "+ele.p())
      }

    }
  }

  reader.close()
}


