import 'holderjs';
import { useState, useEffect } from 'react';
const api=window.api;

const Body= ()=>
{
    const [plotimg,setPlotImg] = useState('holder.js/700x600?text=Plot');
    const [star,setStar] = useState(0);
    const [test,setTest] = useState(1);


    useEffect(() => {
        api.receive("send_plot",(data)=>
        {
            setPlotImg(`data:image/jpeg;base64,${Buffer.from(data).toString("base64")}`);
            setTest(50);
        }); 
    
        return () => {
            
        }

    }, []);

    const getImage = ()=>
    {
        api.send("get_plot",star);
    }

    return(
        
        <div className="container-fluid">
            <div className="row mt-5 justify-content-between border border-dark">
                <div className="col-md-auto">
                    <img src={plotimg} alt="Star plot" width={700} height={600}/>
                </div>

                <div className="col-md-auto">
                    <div>
                        <label htmlFor="star" className="form-label">Star to plot</label>
                        <input type="number" className="form-control mb-3" onChange={(e)=>setStar(e.target.value)}></input>
                        <button type="button" className="btn btn-primary" onClick={()=>{getImage()}}>Select</button>
                        <p>{star}</p>
                        <p>{test}</p>
                    </div>
                </div>
            </div>
        </div>

    )
}

export default Body